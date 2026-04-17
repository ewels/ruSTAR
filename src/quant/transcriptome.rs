//! Transcript-level quantification (`--quantMode TranscriptomeSAM`).
//!
//! Builds a per-transcript exon map from GTF records and projects genome-space
//! alignments onto the transcriptome coordinate system so downstream tools
//! (Salmon, RSEM) can read `Aligned.toTranscriptome.out.bam`.
//!
//! This module mirrors STAR's `Transcriptome` / `alignToTranscript` /
//! `quantAlign` logic (see `source/Transcriptome.cpp` and
//! `source/Transcriptome_quantAlign.cpp` in the upstream repo), with the
//! only substantive divergence that ruSTAR builds the transcript tables on
//! the fly from the input GTF instead of loading persisted
//! `transcriptInfo.tab`/`exonInfo.tab` files.
use std::collections::HashMap;

use crate::align::transcript::{CigarOp, Exon, Transcript};
use crate::error::Error;
use crate::genome::Genome;
use crate::junction::gtf::GtfRecord;

/// Per-transcript exon in absolute genome coordinates (0-based half-open),
/// paired with the cumulative transcript-space length of all preceding exons
/// (STAR's `exLenCum`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TrExon {
    /// Absolute genome start (0-based, inclusive).
    pub genome_start: u64,
    /// Absolute genome end (0-based, exclusive).
    pub genome_end: u64,
    /// Cumulative transcript-space length of all prior exons (0 for the first).
    pub ex_len_cum: u32,
}

/// Transcriptome metadata built from GTF exon records.
///
/// Each transcript is a contiguous set of exons on a single chromosome and
/// strand, ordered by ascending genome start. `tr_order` / `tr_starts_sorted`
/// / `tr_end_max_sorted` support STAR's `binarySearch1a` + running-max
/// early-exit in `quantAlign`.
#[derive(Debug, Clone)]
pub struct TranscriptomeIndex {
    /// Transcript IDs in insertion order (index = transcript_idx).
    pub tr_ids: Vec<String>,
    /// Genome chromosome index per transcript.
    pub tr_chr_idx: Vec<usize>,
    /// Strand per transcript: 1 = forward/`+`, 2 = reverse/`-`.
    pub tr_strand: Vec<u8>,
    /// gene_id string per transcript.
    pub tr_gene_id: Vec<String>,
    /// Absolute genome start of the first exon (0-based).
    pub tr_start: Vec<u64>,
    /// Absolute genome end of the last exon (0-based half-open).
    pub tr_end: Vec<u64>,
    /// Exons per transcript (sorted by genome_start, with `ex_len_cum`).
    pub tr_exons: Vec<Vec<TrExon>>,
    /// Total transcript-space length (sum of exon spans).
    pub tr_length: Vec<u32>,
    /// Permutation of `0..n_transcripts` sorted by `(tr_start, tr_end)`.
    pub tr_order: Vec<usize>,
    /// `tr_start[tr_order[i]]` — for `binary_search`.
    pub tr_starts_sorted: Vec<u64>,
    /// Running max of `tr_end` along `tr_order` (STAR's `trEmax`).
    pub tr_end_max_sorted: Vec<u64>,
}

impl TranscriptomeIndex {
    /// Build from already-parsed GTF exon records.
    ///
    /// Records are grouped by `transcript_id`.  Transcripts with inconsistent
    /// strands / chromosomes are skipped with a warning.  Transcripts on
    /// unknown chromosomes are skipped with a warning.
    pub fn from_gtf_exons(exons: &[GtfRecord], genome: &Genome) -> Result<Self, Error> {
        // Group exons by transcript_id, preserving first-seen insertion order.
        let mut order: Vec<String> = Vec::new();
        let mut groups: HashMap<String, Vec<&GtfRecord>> = HashMap::new();
        for rec in exons {
            let tid = match rec.attributes.get("transcript_id") {
                Some(id) => id.clone(),
                None => continue,
            };
            if !groups.contains_key(&tid) {
                order.push(tid.clone());
            }
            groups.entry(tid).or_default().push(rec);
        }

        let mut tr_ids: Vec<String> = Vec::new();
        let mut tr_chr_idx: Vec<usize> = Vec::new();
        let mut tr_strand: Vec<u8> = Vec::new();
        let mut tr_gene_id: Vec<String> = Vec::new();
        let mut tr_start: Vec<u64> = Vec::new();
        let mut tr_end: Vec<u64> = Vec::new();
        let mut tr_exons_vec: Vec<Vec<TrExon>> = Vec::new();
        let mut tr_length: Vec<u32> = Vec::new();

        for tid in &order {
            let mut recs = groups.remove(tid).unwrap();
            if recs.is_empty() {
                continue;
            }

            // Validate consistent chromosome + strand.
            let first = recs[0];
            let chr_name = &first.seqname;
            let strand_char = first.strand;
            let mut inconsistent = false;
            for r in &recs {
                if &r.seqname != chr_name || r.strand != strand_char {
                    inconsistent = true;
                    break;
                }
            }
            if inconsistent {
                log::warn!(
                    "quantMode TranscriptomeSAM: transcript {} has inconsistent chromosome/strand across exons — skipping",
                    tid
                );
                continue;
            }

            // Map chromosome name to genome index.
            let chr_idx = match genome.chr_name.iter().position(|n| n == chr_name) {
                Some(i) if i < genome.n_chr_real => i,
                _ => {
                    log::warn!(
                        "quantMode TranscriptomeSAM: transcript {} on unknown chromosome {} — skipping",
                        tid,
                        chr_name
                    );
                    continue;
                }
            };
            let chr_offset = genome.chr_start[chr_idx];

            // Sort exons by GTF start ASC.
            recs.sort_by_key(|r| r.start);

            // Build TrExon list in absolute genome coords (0-based half-open).
            let mut tr_exons: Vec<TrExon> = Vec::with_capacity(recs.len());
            let mut ex_len_cum: u32 = 0;
            for r in &recs {
                // GTF is 1-based inclusive; ruSTAR uses 0-based half-open.
                let abs_start = chr_offset + r.start.saturating_sub(1);
                let abs_end = chr_offset + r.end;
                if abs_end <= abs_start {
                    log::warn!(
                        "quantMode TranscriptomeSAM: transcript {} has invalid exon {}-{} — skipping",
                        tid,
                        r.start,
                        r.end
                    );
                    tr_exons.clear();
                    break;
                }
                let len = (abs_end - abs_start) as u32;
                tr_exons.push(TrExon {
                    genome_start: abs_start,
                    genome_end: abs_end,
                    ex_len_cum,
                });
                ex_len_cum = ex_len_cum.saturating_add(len);
            }
            if tr_exons.is_empty() {
                continue;
            }

            let first_ex = tr_exons.first().unwrap();
            let last_ex = tr_exons.last().unwrap();
            let start_abs = first_ex.genome_start;
            let end_abs = last_ex.genome_end;
            let total_len = ex_len_cum;

            let strand_u8 = match strand_char {
                '+' => 1u8,
                '-' => 2u8,
                _ => 1u8, // unknown → treat as forward (STAR default)
            };

            let gene_id = first
                .attributes
                .get("gene_id")
                .cloned()
                .unwrap_or_default();

            tr_ids.push(tid.clone());
            tr_chr_idx.push(chr_idx);
            tr_strand.push(strand_u8);
            tr_gene_id.push(gene_id);
            tr_start.push(start_abs);
            tr_end.push(end_abs);
            tr_exons_vec.push(tr_exons);
            tr_length.push(total_len);
        }

        let n_tr = tr_ids.len();

        // Sorted view by (tr_start, tr_end) for binary-search +
        // running-max early-exit.
        let mut tr_order: Vec<usize> = (0..n_tr).collect();
        tr_order.sort_by(|&a, &b| {
            tr_start[a]
                .cmp(&tr_start[b])
                .then_with(|| tr_end[a].cmp(&tr_end[b]))
        });

        let tr_starts_sorted: Vec<u64> = tr_order.iter().map(|&i| tr_start[i]).collect();
        let mut tr_end_max_sorted: Vec<u64> = Vec::with_capacity(n_tr);
        let mut running_max: u64 = 0;
        for &i in &tr_order {
            running_max = running_max.max(tr_end[i]);
            tr_end_max_sorted.push(running_max);
        }

        Ok(TranscriptomeIndex {
            tr_ids,
            tr_chr_idx,
            tr_strand,
            tr_gene_id,
            tr_start,
            tr_end,
            tr_exons: tr_exons_vec,
            tr_length,
            tr_order,
            tr_starts_sorted,
            tr_end_max_sorted,
        })
    }

    /// Number of transcripts indexed.
    pub fn n_transcripts(&self) -> usize {
        self.tr_ids.len()
    }
}

// ---------------------------------------------------------------------------
// Projection (subtask 2) — stub; implemented in a later commit.
// ---------------------------------------------------------------------------

/// Project a genome-space `Transcript` onto every transcript whose coordinates
/// fully contain the alignment.  Returns transcript-space alignments ready for
/// BAM emission.
///
/// Implementation lives in subtask 2.
pub fn align_to_transcripts(
    _align: &Transcript,
    _idx: &TranscriptomeIndex,
    _lread: u32,
) -> Vec<Transcript> {
    // Filled in during subtask 2.
    Vec::new()
}

#[allow(dead_code)]
#[doc(hidden)]
fn _touch_unused_imports(_cigar: &CigarOp, _exon: &Exon) {}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap as StdHashMap;

    fn make_genome() -> Genome {
        Genome {
            sequence: vec![0u8; 3000],
            n_genome: 3000,
            n_chr_real: 2,
            chr_start: vec![0, 1000, 3000],
            chr_length: vec![1000, 2000],
            chr_name: vec!["chr1".to_string(), "chr2".to_string()],
        }
    }

    fn make_exon(
        seqname: &str,
        start: u64,
        end: u64,
        strand: char,
        gene_id: &str,
        transcript_id: &str,
    ) -> GtfRecord {
        let mut attrs: StdHashMap<String, String> = StdHashMap::new();
        attrs.insert("gene_id".to_string(), gene_id.to_string());
        attrs.insert("transcript_id".to_string(), transcript_id.to_string());
        GtfRecord {
            seqname: seqname.to_string(),
            feature: "exon".to_string(),
            start,
            end,
            strand,
            attributes: attrs,
        }
    }

    #[test]
    fn single_exon_transcript_metadata() {
        let genome = make_genome();
        let exons = vec![make_exon("chr1", 101, 200, '+', "G1", "T1")];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        assert_eq!(idx.n_transcripts(), 1);
        assert_eq!(idx.tr_ids[0], "T1");
        assert_eq!(idx.tr_gene_id[0], "G1");
        assert_eq!(idx.tr_strand[0], 1);
        assert_eq!(idx.tr_chr_idx[0], 0);
        // GTF 1-based 101 → absolute 0-based 100; end 200 → absolute 200 (exclusive)
        assert_eq!(idx.tr_start[0], 100);
        assert_eq!(idx.tr_end[0], 200);
        assert_eq!(idx.tr_length[0], 100);
        assert_eq!(idx.tr_exons[0].len(), 1);
        assert_eq!(idx.tr_exons[0][0].ex_len_cum, 0);
    }

    #[test]
    fn multi_exon_transcript_ex_len_cum() {
        let genome = make_genome();
        let exons = vec![
            make_exon("chr1", 101, 200, '+', "G1", "T1"),   // len 100
            make_exon("chr1", 301, 400, '+', "G1", "T1"),   // len 100
            make_exon("chr1", 501, 650, '+', "G1", "T1"),   // len 150
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        assert_eq!(idx.n_transcripts(), 1);
        assert_eq!(idx.tr_length[0], 350);

        let ex = &idx.tr_exons[0];
        assert_eq!(ex[0].ex_len_cum, 0);
        assert_eq!(ex[1].ex_len_cum, 100);
        assert_eq!(ex[2].ex_len_cum, 200);

        // Ex_len_cum must be monotonically non-decreasing.
        for w in ex.windows(2) {
            assert!(w[0].ex_len_cum <= w[1].ex_len_cum);
        }

        // tr_start / tr_end from first/last absolute exon.
        assert_eq!(idx.tr_start[0], 100);
        assert_eq!(idx.tr_end[0], 650);
    }

    #[test]
    fn reverse_strand_transcript() {
        let genome = make_genome();
        let exons = vec![
            make_exon("chr1", 101, 200, '-', "G2", "T2"),
            make_exon("chr1", 301, 400, '-', "G2", "T2"),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        assert_eq!(idx.n_transcripts(), 1);
        assert_eq!(idx.tr_strand[0], 2);
        assert_eq!(idx.tr_length[0], 200);
    }

    #[test]
    fn two_transcripts_same_gene() {
        let genome = make_genome();
        let exons = vec![
            make_exon("chr1", 101, 200, '+', "G1", "T1"),
            make_exon("chr1", 301, 400, '+', "G1", "T1"),
            make_exon("chr1", 101, 200, '+', "G1", "T2"),
            make_exon("chr1", 301, 500, '+', "G1", "T2"),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        assert_eq!(idx.n_transcripts(), 2);
        assert_eq!(idx.tr_gene_id[0], "G1");
        assert_eq!(idx.tr_gene_id[1], "G1");
        // Different lengths: T1 = 100+100 = 200, T2 = 100+200 = 300
        assert_eq!(idx.tr_length[0], 200);
        assert_eq!(idx.tr_length[1], 300);
    }

    #[test]
    fn unknown_chromosome_skipped() {
        let genome = make_genome();
        let exons = vec![make_exon("chrX", 101, 200, '+', "G1", "T1")];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        assert_eq!(idx.n_transcripts(), 0);
    }

    #[test]
    fn inconsistent_strand_skipped() {
        let genome = make_genome();
        let exons = vec![
            make_exon("chr1", 101, 200, '+', "G1", "T1"),
            make_exon("chr1", 301, 400, '-', "G1", "T1"),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        assert_eq!(idx.n_transcripts(), 0);
    }

    #[test]
    fn tr_end_max_sorted_is_running_max() {
        let genome = make_genome();
        let exons = vec![
            // T1 chr1 starts at 100, ends at 200
            make_exon("chr1", 101, 200, '+', "G1", "T1"),
            // T2 chr1 starts at 150, ends at 500 (encloses T1)
            make_exon("chr1", 151, 500, '+', "G1", "T2"),
            // T3 chr1 starts at 200, ends at 300 (nested)
            make_exon("chr1", 201, 300, '+', "G1", "T3"),
        ];
        let idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        assert_eq!(idx.n_transcripts(), 3);

        // tr_order must be sorted by (start, end)
        let sorted_starts: Vec<u64> = idx.tr_starts_sorted.clone();
        let mut check = sorted_starts.clone();
        check.sort();
        assert_eq!(sorted_starts, check);

        // tr_end_max_sorted must be monotonically non-decreasing
        for w in idx.tr_end_max_sorted.windows(2) {
            assert!(w[0] <= w[1]);
        }
        // Last entry equals overall maximum tr_end
        let overall_max = *idx.tr_end.iter().max().unwrap();
        assert_eq!(*idx.tr_end_max_sorted.last().unwrap(), overall_max);
    }
}
