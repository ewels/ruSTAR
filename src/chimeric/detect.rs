// Chimeric alignment detection algorithms

use crate::align::SeedCluster;
use crate::align::stitch::stitch_seeds;
use crate::align::transcript::Transcript;
use crate::chimeric::score::{calculate_repeat_length, classify_junction_type};
use crate::chimeric::segment::{ChimericAlignment, ChimericSegment};
use crate::error::Error;
use crate::index::GenomeIndex;
use crate::params::Parameters;

/// Chimeric alignment detector
pub struct ChimericDetector<'a> {
    params: &'a Parameters,
}

impl<'a> ChimericDetector<'a> {
    /// Create a new chimeric detector
    pub fn new(params: &'a Parameters) -> Self {
        Self { params }
    }

    /// Detect chimeric alignments from soft-clipped transcripts (Tier 1)
    ///
    /// Triggers when:
    /// - Transcript has >20% soft-clipped bases
    /// - Soft-clipped segment >= chimSegmentMin
    pub fn detect_from_soft_clips(
        &self,
        transcript: &Transcript,
        read_seq: &[u8],
        _read_name: &str,
        _index: &GenomeIndex,
    ) -> Result<Option<ChimericAlignment>, Error> {
        let (left_clip, right_clip) = transcript.count_soft_clips();
        let total_clip = left_clip + right_clip;
        let read_len = read_seq.len() as u32;

        // Check if soft-clipping exceeds threshold
        if total_clip as f64 / read_len as f64 <= 0.20 {
            return Ok(None);
        }

        // Check if clipped segment meets minimum length
        let min_seg = self.params.chim_segment_min;
        if left_clip < min_seg && right_clip < min_seg {
            return Ok(None);
        }

        // For now, we mark this as a potential chimera but don't re-map the clipped region
        // (Tier 3 re-mapping is deferred to Phase 12.2)
        // Instead, we just log that a chimeric candidate was detected but not processed
        Ok(None)
    }

    /// Detect chimeric alignments from multi-cluster seeds (Tier 2)
    ///
    /// Triggers when:
    /// - Seeds cluster on different chromosomes
    /// - Seeds cluster on different strands (same chromosome)
    /// - Seeds cluster with large genomic distance (>1Mb, same chr/strand)
    pub fn detect_from_multi_clusters(
        &self,
        clusters: &[SeedCluster],
        read_seq: &[u8],
        read_name: &str,
        index: &GenomeIndex,
    ) -> Result<Vec<ChimericAlignment>, Error> {
        let mut chimeras = Vec::new();

        // Find cluster pairs with chimeric signatures
        for i in 0..clusters.len() {
            for j in (i + 1)..clusters.len() {
                if self.is_chimeric_signature(&clusters[i], &clusters[j]) {
                    // Try to build chimeric alignment from these clusters
                    if let Some(chim) = self.build_chimeric_from_clusters(
                        &clusters[i],
                        &clusters[j],
                        read_seq,
                        read_name,
                        index,
                    )? {
                        chimeras.push(chim);
                    }
                }
            }
        }

        Ok(chimeras)
    }

    /// Check if two clusters represent a chimeric signature
    fn is_chimeric_signature(&self, c1: &SeedCluster, c2: &SeedCluster) -> bool {
        // Different chromosomes
        if c1.chr_idx != c2.chr_idx {
            return true;
        }

        // Different strands (same chromosome)
        if c1.is_reverse != c2.is_reverse {
            return true;
        }

        // Large genomic distance (same chr/strand)
        let distance = genomic_distance(c1, c2);
        if distance > 1_000_000 {
            return true;
        }

        false
    }

    /// Build chimeric alignment from two clusters
    fn build_chimeric_from_clusters(
        &self,
        cluster1: &SeedCluster,
        cluster2: &SeedCluster,
        read_seq: &[u8],
        read_name: &str,
        index: &GenomeIndex,
    ) -> Result<Option<ChimericAlignment>, Error> {
        if cluster1.alignments.is_empty() || cluster2.alignments.is_empty() {
            return Ok(None);
        }

        // Stitch each cluster independently using existing stitch_seeds
        use crate::align::score::AlignmentScorer;
        let scorer = AlignmentScorer::from_params(self.params);

        let transcripts1 = stitch_seeds(cluster1, read_seq, index, &scorer)?;
        let transcripts2 = stitch_seeds(cluster2, read_seq, index, &scorer)?;

        if transcripts1.is_empty() || transcripts2.is_empty() {
            return Ok(None);
        }

        // Take best transcript from each cluster
        let t1 = &transcripts1[0];
        let t2 = &transcripts2[0];

        // Check that both transcripts have exons
        if t1.exons.is_empty() || t2.exons.is_empty() {
            return Ok(None);
        }

        // Determine donor/acceptor based on read position
        let (donor_t, acceptor_t) = if t1.exons[0].read_start < t2.exons[0].read_start {
            (t1, t2)
        } else {
            (t2, t1)
        };

        // Convert transcripts to chimeric segments
        let donor = transcript_to_segment(donor_t)?;
        let acceptor = transcript_to_segment(acceptor_t)?;

        // Check minimum segment lengths
        if !donor.meets_min_length(self.params.chim_segment_min)
            || !acceptor.meets_min_length(self.params.chim_segment_min)
        {
            return Ok(None);
        }

        // Classify junction type
        let junction_type = classify_junction_type(
            &index.genome,
            donor.chr_idx,
            donor.genome_end,
            donor.is_reverse,
            acceptor.chr_idx,
            acceptor.genome_start,
            acceptor.is_reverse,
        );

        // Calculate repeat lengths
        let (repeat_len_donor, repeat_len_acceptor) = calculate_repeat_length(
            &index.genome,
            donor.chr_idx,
            donor.genome_end,
            acceptor.chr_idx,
            acceptor.genome_start,
            20, // max check distance
        );

        // Create chimeric alignment
        let chim = ChimericAlignment::new(
            donor,
            acceptor,
            junction_type,
            repeat_len_donor,
            repeat_len_acceptor,
            read_seq.to_vec(),
            read_name.to_string(),
        );

        Ok(Some(chim))
    }
}

/// Calculate genomic distance between two clusters
fn genomic_distance(c1: &SeedCluster, c2: &SeedCluster) -> u64 {
    if c1.chr_idx != c2.chr_idx {
        return u64::MAX;
    }

    if c1.genome_end < c2.genome_start {
        c2.genome_start - c1.genome_end
    } else {
        c1.genome_start.saturating_sub(c2.genome_end)
    }
}

/// Detect inter-mate chimeric alignment from two single-mate transcripts.
///
/// Fires when mate1 and mate2 map to different chromosomes, opposite-orientation
/// same-chromosome positions, or positions too far apart to be a normal PE pair.
/// This is the primary PE-specific chimeric case (gene-fusion detection).
pub fn detect_inter_mate_chimeric(
    t1: &Transcript,
    t2: &Transcript,
    mate1_seq: &[u8],
    read_name: &str,
    params: &Parameters,
    index: &GenomeIndex,
) -> Option<ChimericAlignment> {
    // Only fire if the pair is discordant (different chr, same strand (both FW or both RC =
    // not FR orientation), or too far apart).
    let is_inter_chr = t1.chr_idx != t2.chr_idx;
    // FR pair expects t1.is_reverse=false (mate1 FW) and t2.is_reverse=true (mate2 RC).
    // Chimeric if both same strand.
    let same_strand = t1.is_reverse == t2.is_reverse;
    let too_far = if t1.chr_idx == t2.chr_idx {
        let left_end = t1.genome_end.min(t2.genome_end);
        let right_start = t1.genome_start.max(t2.genome_start);
        right_start > left_end && right_start - left_end > 1_000_000
    } else {
        false
    };

    if !is_inter_chr && !same_strand && !too_far {
        return None;
    }

    if t1.exons.is_empty() || t2.exons.is_empty() {
        return None;
    }

    // Convert transcripts to chimeric segments
    let donor = transcript_to_segment(t1).ok()?;
    let acceptor = transcript_to_segment(t2).ok()?;

    if !donor.meets_min_length(params.chim_segment_min)
        || !acceptor.meets_min_length(params.chim_segment_min)
    {
        return None;
    }

    // Junction type: non-canonical (0) for inter-chromosomal; try motif for same-chr
    let junction_type = if is_inter_chr {
        0
    } else {
        classify_junction_type(
            &index.genome,
            donor.chr_idx,
            donor.genome_end,
            donor.is_reverse,
            acceptor.chr_idx,
            acceptor.genome_start,
            acceptor.is_reverse,
        )
    };

    let (repeat_len_donor, repeat_len_acceptor) = calculate_repeat_length(
        &index.genome,
        donor.chr_idx,
        donor.genome_end,
        acceptor.chr_idx,
        acceptor.genome_start,
        20,
    );

    let chim = ChimericAlignment::new(
        donor,
        acceptor,
        junction_type,
        repeat_len_donor,
        repeat_len_acceptor,
        mate1_seq.to_vec(),
        read_name.to_string(),
    );

    Some(chim)
}

/// Convert a transcript to a chimeric segment
pub(crate) fn transcript_to_segment(transcript: &Transcript) -> Result<ChimericSegment, Error> {
    if transcript.exons.is_empty() {
        return Err(Error::Alignment(
            "Cannot convert empty transcript to segment".to_string(),
        ));
    }

    // Get overall bounds
    let read_start = transcript.exons[0].read_start;
    let read_end = transcript.exons.last().unwrap().read_end;

    Ok(ChimericSegment {
        chr_idx: transcript.chr_idx,
        genome_start: transcript.genome_start,
        genome_end: transcript.genome_end,
        is_reverse: transcript.is_reverse,
        read_start,
        read_end,
        cigar: transcript.cigar.clone(),
        score: transcript.score,
        n_mismatch: transcript.n_mismatch,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::WindowAlignment;
    use crate::align::transcript::{CigarOp, Exon, Transcript};
    use crate::genome::Genome;
    use crate::index::GenomeIndex;
    use crate::index::packed_array::PackedArray;
    use crate::index::sa_index::SaIndex;
    use crate::index::suffix_array::SuffixArray;
    use crate::junction::SpliceJunctionDb;
    use clap::Parser;

    /// Minimal two-chromosome genome for chimeric tests.
    /// Each chromosome has 200 bases of A (=0), padded to 256-byte bins.
    fn make_test_genome() -> Genome {
        let chr_len = 200u64;
        let chr_pad = 256u64;
        let n_genome = chr_pad * 2;
        let sequence = vec![0u8; 2 * n_genome as usize];
        Genome {
            sequence,
            n_genome,
            n_chr_real: 2,
            chr_name: vec!["chr0".to_string(), "chr1".to_string()],
            chr_length: vec![chr_len, chr_len],
            chr_start: vec![0, chr_pad, n_genome],
        }
    }

    fn make_test_index() -> GenomeIndex {
        let genome = make_test_genome();
        let gstrand_bit = 32u32;
        GenomeIndex {
            genome,
            suffix_array: SuffixArray {
                data: PackedArray::new(33, 0),
                gstrand_bit,
                gstrand_mask: (1u64 << gstrand_bit) - 1,
            },
            sa_index: SaIndex {
                nbases: 0,
                genome_sa_index_start: vec![0],
                data: PackedArray::new(35, 0),
                word_length: 35,
                gstrand_bit,
            },
            junction_db: SpliceJunctionDb::empty(),
            transcriptome: None,
            prepared_junctions: Vec::new(),
        }
    }

    fn make_transcript(
        chr_idx: usize,
        genome_start: u64,
        genome_end: u64,
        is_reverse: bool,
    ) -> Transcript {
        let read_len = (genome_end - genome_start) as usize;
        Transcript {
            chr_idx,
            genome_start,
            genome_end,
            is_reverse,
            exons: vec![Exon {
                genome_start,
                genome_end,
                read_start: 0,
                read_end: read_len,
                i_frag: 0,
            }],
            cigar: vec![CigarOp::Match(read_len as u32)],
            score: read_len as i32,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![0u8; read_len],
        }
    }

    /// Helper to create a minimal SeedCluster for chimeric detection tests
    fn make_test_cluster(
        chr_idx: usize,
        genome_start: u64,
        genome_end: u64,
        is_reverse: bool,
    ) -> SeedCluster {
        SeedCluster {
            alignments: vec![WindowAlignment {
                seed_idx: 0,
                read_pos: 0,
                length: (genome_end - genome_start) as usize,
                genome_pos: genome_start,
                sa_pos: genome_start,
                n_rep: 1,
                is_anchor: true,
                mate_id: 2,
                pre_ext_score: (genome_end - genome_start) as i32,
            }],
            chr_idx,
            genome_start,
            genome_end,
            is_reverse,
            anchor_idx: 0,
            anchor_bin: 0,
        }
    }

    #[test]
    fn test_genomic_distance_same_chr() {
        let c1 = make_test_cluster(0, 1000, 1100, false);
        let c2 = make_test_cluster(0, 1200, 1300, false);

        assert_eq!(genomic_distance(&c1, &c2), 100);
        assert_eq!(genomic_distance(&c2, &c1), 100);
    }

    #[test]
    fn test_genomic_distance_overlapping() {
        let c1 = make_test_cluster(0, 1000, 1200, false);
        let c2 = make_test_cluster(0, 1100, 1300, false);

        assert_eq!(genomic_distance(&c1, &c2), 0);
    }

    #[test]
    fn test_genomic_distance_different_chr() {
        let c1 = make_test_cluster(0, 1000, 1100, false);
        let c2 = make_test_cluster(1, 1000, 1100, false);

        assert_eq!(genomic_distance(&c1, &c2), u64::MAX);
    }

    #[test]
    fn test_is_chimeric_signature_different_chr() {
        let params = Parameters::try_parse_from(vec!["ruSTAR"]).unwrap();
        let detector = ChimericDetector::new(&params);

        let c1 = make_test_cluster(0, 1000, 1100, false);
        let c2 = make_test_cluster(1, 1000, 1100, false);

        assert!(detector.is_chimeric_signature(&c1, &c2));
    }

    #[test]
    fn test_is_chimeric_signature_strand_break() {
        let params = Parameters::try_parse_from(vec!["ruSTAR"]).unwrap();
        let detector = ChimericDetector::new(&params);

        let c1 = make_test_cluster(0, 1000, 1100, false);
        let c2 = make_test_cluster(0, 1200, 1300, true);

        assert!(detector.is_chimeric_signature(&c1, &c2));
    }

    #[test]
    fn test_is_chimeric_signature_large_distance() {
        let params = Parameters::try_parse_from(vec!["ruSTAR"]).unwrap();
        let detector = ChimericDetector::new(&params);

        let c1 = make_test_cluster(0, 1000, 1100, false);
        let c2 = make_test_cluster(0, 2_000_000, 2_000_100, false);

        assert!(detector.is_chimeric_signature(&c1, &c2));
    }

    #[test]
    fn test_is_chimeric_signature_close_same_strand() {
        let params = Parameters::try_parse_from(vec!["ruSTAR"]).unwrap();
        let detector = ChimericDetector::new(&params);

        let c1 = make_test_cluster(0, 1000, 1100, false);
        let c2 = make_test_cluster(0, 1200, 1300, false);

        assert!(!detector.is_chimeric_signature(&c1, &c2));
    }

    // --- transcript_to_segment tests ---

    #[test]
    fn test_transcript_to_segment_basic() {
        let t = make_transcript(0, 1000, 1100, false);
        let seg = transcript_to_segment(&t).unwrap();

        assert_eq!(seg.chr_idx, 0);
        assert_eq!(seg.genome_start, 1000);
        assert_eq!(seg.genome_end, 1100);
        assert!(!seg.is_reverse);
        assert_eq!(seg.read_start, 0);
        assert_eq!(seg.read_end, 100);
        assert_eq!(seg.score, 100);
    }

    #[test]
    fn test_transcript_to_segment_empty_returns_error() {
        let t = Transcript {
            chr_idx: 0,
            genome_start: 0,
            genome_end: 0,
            is_reverse: false,
            exons: vec![],
            cigar: vec![],
            score: 0,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![],
        };
        assert!(transcript_to_segment(&t).is_err());
    }

    // --- detect_inter_mate_chimeric tests ---

    #[test]
    fn test_inter_mate_chimeric_concordant_returns_none() {
        // Normal FR pair on the same chromosome, close together → not chimeric
        let params = Parameters::try_parse_from(vec!["ruSTAR", "--chimSegmentMin", "10"]).unwrap();
        let index = make_test_index();

        let t1 = make_transcript(0, 10, 60, false); // mate1 forward
        let t2 = make_transcript(0, 80, 130, true); // mate2 reverse, same chr, close
        let read_seq = vec![0u8; 50];

        let result = detect_inter_mate_chimeric(&t1, &t2, &read_seq, "read1", &params, &index);
        assert!(result.is_none());
    }

    #[test]
    fn test_inter_mate_chimeric_different_chromosomes() {
        let params = Parameters::try_parse_from(vec!["ruSTAR", "--chimSegmentMin", "10"]).unwrap();
        let index = make_test_index();

        let t1 = make_transcript(0, 10, 60, false); // mate1 chr0
        let t2 = make_transcript(1, 10, 60, true); // mate2 chr1
        let read_seq = vec![0u8; 50];

        let result = detect_inter_mate_chimeric(&t1, &t2, &read_seq, "read1", &params, &index);
        assert!(result.is_some());
        let chim = result.unwrap();
        // Donor is the mate with earlier read_start (both 0 here; donor is t1 by read_start tie)
        assert_ne!(chim.donor.chr_idx, chim.acceptor.chr_idx);
    }

    #[test]
    fn test_inter_mate_chimeric_same_strand() {
        // Both mates forward on the same chromosome → chimeric (strand break)
        let params = Parameters::try_parse_from(vec!["ruSTAR", "--chimSegmentMin", "10"]).unwrap();
        let index = make_test_index();

        let t1 = make_transcript(0, 10, 60, false); // mate1 forward
        let t2 = make_transcript(0, 80, 130, false); // mate2 also forward (abnormal)
        let read_seq = vec![0u8; 50];

        let result = detect_inter_mate_chimeric(&t1, &t2, &read_seq, "read1", &params, &index);
        assert!(result.is_some());
    }

    #[test]
    fn test_inter_mate_chimeric_too_far() {
        // Opposite-strand pair but >1Mb apart → chimeric
        let params = Parameters::try_parse_from(vec!["ruSTAR", "--chimSegmentMin", "10"]).unwrap();
        let index = make_test_index();

        // Use large positions — out-of-bounds for sequence but score.rs guards handle this
        let t1 = make_transcript(0, 10, 60, false);
        let t2 = make_transcript(0, 2_000_000, 2_000_050, true);
        let read_seq = vec![0u8; 50];

        let result = detect_inter_mate_chimeric(&t1, &t2, &read_seq, "read1", &params, &index);
        assert!(result.is_some());
    }

    #[test]
    fn test_inter_mate_chimeric_segment_too_short() {
        // chimSegmentMin=100 but segments are only 20bp → None
        let params = Parameters::try_parse_from(vec!["ruSTAR", "--chimSegmentMin", "100"]).unwrap();
        let index = make_test_index();

        let t1 = make_transcript(0, 10, 30, false);
        let t2 = make_transcript(1, 10, 30, true);
        let read_seq = vec![0u8; 20];

        let result = detect_inter_mate_chimeric(&t1, &t2, &read_seq, "read1", &params, &index);
        assert!(result.is_none());
    }

    #[test]
    fn test_inter_mate_chimeric_empty_exons_returns_none() {
        let params = Parameters::try_parse_from(vec!["ruSTAR", "--chimSegmentMin", "10"]).unwrap();
        let index = make_test_index();
        let read_seq = vec![0u8; 50];

        let t1 = make_transcript(0, 10, 60, false);
        let mut t2 = make_transcript(1, 10, 60, true);
        t2.exons.clear();

        let result = detect_inter_mate_chimeric(&t1, &t2, &read_seq, "read1", &params, &index);
        assert!(result.is_none());
    }
}
