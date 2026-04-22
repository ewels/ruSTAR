//! STAR-faithful splice junction insertion into the genome index.
//!
//! Ports `source/sjdbPrepare.cpp` + `source/sjdbBuildIndex.cpp` from STAR.
//! At `genomeGenerate` time, STAR extracts the flanking `sjdbOverhang`
//! bases on each side of every GTF-derived splice junction, concatenates
//! them into a `Gsj` buffer, appends that buffer to the `Genome` binary,
//! and extends the suffix array to index the new bases. This module
//! provides the same machinery for ruSTAR so that the generated
//! `Genome` / `SA` / `SAindex` / `sjdbInfo.txt` / `sjdbList.out.tab` files
//! match STAR's byte-for-byte.
//!
//! The module is orchestrated from `index::GenomeIndex::build` after the
//! base suffix array has been built.

use crate::align::score::detect_splice_motif;
use crate::genome::Genome;
use crate::junction::encode_motif;

/// Compute STAR's `(sjdbShiftLeft, sjdbShiftRight)` for an intron whose
/// 0-based donor/acceptor positions are `s` and `e`.
///
/// STAR defines these as the number of bases the intron boundary can shift
/// left / right while preserving the donor/acceptor base identity
/// (`sjdbPrepare.cpp:52-73`) — i.e. the repeat length across the junction.
/// Intended to land the junction at its left-most canonical position so
/// identical splice events produce identical SA indices regardless of
/// which exon pair they came from.
///
/// Stops at genome bounds, on any N-base (code ≥ 4), or at the 255 cap.
pub fn compute_shifts(genome: &Genome, s: u64, e: u64, n_genome_real: u64) -> (u8, u8) {
    let forward = &genome.sequence[..n_genome_real as usize];
    let si = s as usize;
    let ei = e as usize;

    let mut jj_l: u8 = 0;
    while jj_l < 255 && (jj_l as usize) < si && ei >= jj_l as usize {
        let a = forward[si - 1 - jj_l as usize];
        let b = forward[ei - jj_l as usize];
        if a != b || a >= 4 {
            break;
        }
        jj_l += 1;
    }

    let mut jj_r: u8 = 0;
    while jj_r < 255 && (ei + 1 + jj_r as usize) < forward.len() {
        let a = forward[si + jj_r as usize];
        let b = forward[ei + 1 + jj_r as usize];
        if a != b || a >= 4 {
            break;
        }
        jj_r += 1;
    }

    (jj_l, jj_r)
}

/// A single junction with all metadata STAR needs to write the sjdb
/// files and extend the genome. Positions are 0-based absolute genome
/// coordinates; `start_pos` / `end_pos` are the FIRST and LAST bases of
/// the intron (inclusive), already shifted left by `shift_left` so they
/// sit at the canonical (left-most) representation of the motif.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PreparedJunction {
    /// Chromosome index the junction belongs to.
    pub chr_idx: usize,
    /// Shift-adjusted 0-based genome position of the first intron base.
    pub start_pos: u64,
    /// Shift-adjusted 0-based genome position of the last intron base.
    pub end_pos: u64,
    /// STAR motif code (0 = non-canonical, 1-6 = canonical variants).
    pub motif: u8,
    /// Repeat length to the left of the (pre-shift) donor.
    pub shift_left: u8,
    /// Repeat length to the right of the (pre-shift) acceptor.
    pub shift_right: u8,
    /// STAR strand code (0 = unknown/dot, 1 = +, 2 = -).
    pub strand: u8,
}

impl PreparedJunction {
    /// Value STAR writes into `mapGen.sjdbStart` after its post-dedup
    /// sort: ORIGINAL pre-shift start for canonical motifs, shifted
    /// start for non-canonical (matches `sjdbPrepare.cpp:127,174`).
    pub fn stored_start(&self) -> u64 {
        if self.motif == 0 {
            self.start_pos
        } else {
            self.start_pos + self.shift_left as u64
        }
    }

    /// Companion to `stored_start` — the `mapGen.sjdbEnd` value.
    pub fn stored_end(&self) -> u64 {
        if self.motif == 0 {
            self.end_pos
        } else {
            self.end_pos + self.shift_left as u64
        }
    }

    /// Original (pre-shift) 0-based position of the first intron base.
    /// Used for Gsj flanking-sequence extraction (STAR uses this regardless
    /// of motif).
    pub fn original_start(&self) -> u64 {
        self.start_pos + self.shift_left as u64
    }

    /// Companion to `original_start`.
    pub fn original_end(&self) -> u64 {
        self.end_pos + self.shift_left as u64
    }
}

/// Convert a splice-junction database entry to a fully-prepared entry
/// carrying motif, shifts, strand, and shift-adjusted coordinates.
///
/// `db_strand` is STAR's 0/1/2 (unknown/+/-); when it's 0 and the motif
/// is canonical, STAR derives the strand from the motif via
/// `2 - motif % 2` (see `sjdbPrepare.cpp:184-188`).
pub fn prepare_junction(
    chr_idx: usize,
    intron_start: u64,
    intron_end: u64,
    db_strand: u8,
    genome: &Genome,
    n_genome_real: u64,
) -> PreparedJunction {
    let intron_len = (intron_end - intron_start + 1) as u32;
    let motif = encode_motif(detect_splice_motif(intron_start, intron_len, genome));
    let (shift_left, shift_right) = compute_shifts(genome, intron_start, intron_end, n_genome_real);

    // sjdbPrepare.cpp:71-72 — land the junction at its left-most canonical
    // representation so identical splice events produce identical indices.
    let shifted_start = intron_start - shift_left as u64;
    let shifted_end = intron_end - shift_left as u64;

    let strand = match db_strand {
        1 | 2 => db_strand,
        _ if motif == 0 => 0,
        _ => 2 - (motif % 2), // 1/3/5 → 1 (+), 2/4/6 → 2 (-)
    };

    PreparedJunction {
        chr_idx,
        start_pos: shifted_start,
        end_pos: shifted_end,
        motif,
        shift_left,
        shift_right,
        strand,
    }
}

/// Sort a prepared junction list into STAR's post-dedup order and apply
/// the cross-strand deduplication that STAR does after its second sort
/// (`sjdbPrepare.cpp:141-192`).
///
/// STAR's first-pass (intra-strand) dedup collapses duplicate sjdb
/// entries from the same source at the same `(start, end, strand)`.
/// ruSTAR's `SpliceJunctionDb` already deduplicates on that key at the
/// HashMap level, so those first-pass branches never trigger here; the
/// second-pass cross-strand collision dedup does.
///
/// Dedup rules when two surviving junctions share `(stored_start,
/// stored_end)` but have different strand assignments:
///
/// - Undefined strand vs defined strand → keep the defined-strand one.
/// - Both non-canonical → collapse to a single entry with strand = 0
///   (undefined).
/// - One canonical + one not → keep the canonical one.
/// - Both canonical but on correct vs wrong strand relative to motif —
///   keep the one whose strand matches `2 - motif % 2`.
pub fn sort_and_dedup(mut junctions: Vec<PreparedJunction>) -> Vec<PreparedJunction> {
    junctions.sort_by(|a, b| {
        a.stored_start()
            .cmp(&b.stored_start())
            .then_with(|| a.stored_end().cmp(&b.stored_end()))
    });

    let mut out: Vec<PreparedJunction> = Vec::with_capacity(junctions.len());
    for j in junctions {
        match out.last() {
            Some(last)
                if last.stored_start() == j.stored_start()
                    && last.stored_end() == j.stored_end() =>
            {
                if let Some(winner) = merge_cross_strand(last, &j) {
                    *out.last_mut().unwrap() = winner;
                }
                // else: keep `last` unchanged.
            }
            _ => out.push(j),
        }
    }
    out
}

/// Decide what `(stored_start, stored_end)` duplicate to keep.
/// Returns `Some(new)` to replace the stored entry, `None` to keep it.
/// For the "both non-canonical on opposite strands" case we keep the
/// existing entry but force its strand to 0 in-place (handled by the
/// caller via a special-case — represented here as returning a cloned
/// `old` with strand=0).
fn merge_cross_strand(
    old: &PreparedJunction,
    new: &PreparedJunction,
) -> Option<PreparedJunction> {
    // Strand 0 = undefined.
    if old.strand > 0 && new.strand == 0 {
        return None; // keep old
    }
    if old.strand == 0 && new.strand > 0 {
        return Some(new.clone()); // replace
    }
    // Both non-canonical → collapse to undefined strand on the old one.
    if old.motif == 0 && new.motif == 0 {
        let mut merged = old.clone();
        merged.strand = 0;
        return Some(merged);
    }
    // One canonical, one not: prefer canonical.
    if old.motif > 0 && new.motif == 0 {
        return None;
    }
    if old.motif == 0 && new.motif > 0 {
        return Some(new.clone());
    }
    // Both canonical with defined strands. Keep the one on the correct
    // strand for its motif (2 - motif % 2). If the old one is on the
    // correct strand, skip the new one; otherwise replace.
    let old_expected = 2 - (old.motif % 2);
    if old.strand == old_expected {
        None
    } else {
        Some(new.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Build a single-chromosome Genome from a forward-strand byte slice.
    // Callers who only exercise compute_shifts / prepare_junction don't
    // need a valid reverse complement; we fill the RC half with padding.
    fn make_test_genome(forward: Vec<u8>) -> Genome {
        let n = forward.len() as u64;
        let mut seq = vec![5u8; (n * 2) as usize];
        seq[..forward.len()].copy_from_slice(&forward);
        Genome {
            sequence: seq,
            n_genome: n,
            n_chr_real: 1,
            chr_name: vec!["chr1".to_string()],
            chr_length: vec![n],
            chr_start: vec![0, n],
        }
    }

    #[test]
    fn shifts_no_repeat() {
        let mut f = vec![4u8; 300];
        f[100] = 0;
        f[101] = 1;
        f[102] = 0;
        f[103] = 1;
        f[104] = 0;
        f[105] = 2; // intron start
        f[106] = 3;
        f[203] = 0;
        f[204] = 2; // intron end
        f[205] = 3;
        f[206] = 2;
        f[207] = 3;
        f[208] = 2;
        f[209] = 3;
        let g = make_test_genome(f);
        let (l, r) = compute_shifts(&g, 105, 204, g.n_genome);
        assert_eq!(l, 0);
        assert_eq!(r, 0);
    }

    #[test]
    fn shifts_with_repeat_on_left() {
        let mut f = vec![4u8; 300];
        f[103] = 0;
        f[104] = 1;
        f[105] = 2;
        f[106] = 3;
        f[203] = 0;
        f[204] = 1;
        f[205] = 3;
        let g = make_test_genome(f);
        let (l, _r) = compute_shifts(&g, 105, 204, g.n_genome);
        assert_eq!(l, 2);
    }

    #[test]
    fn shifts_with_repeat_on_right() {
        let mut f = vec![4u8; 300];
        f[105] = 2;
        f[106] = 3;
        f[107] = 0;
        f[204] = 1;
        f[205] = 2;
        f[206] = 3;
        f[207] = 1;
        let g = make_test_genome(f);
        let (_l, r) = compute_shifts(&g, 105, 204, g.n_genome);
        assert_eq!(r, 2);
    }

    #[test]
    fn shifts_cap_at_255() {
        let g = make_test_genome(vec![0u8; 2000]);
        let (l, r) = compute_shifts(&g, 500, 1000, g.n_genome);
        assert_eq!(l, 255);
        assert_eq!(r, 255);
    }

    #[test]
    fn shifts_stop_at_n_base() {
        let mut f = vec![0u8; 300];
        f[100] = 4;
        let g = make_test_genome(f);
        let (l, _) = compute_shifts(&g, 150, 249, g.n_genome);
        assert_eq!(l, 49);
    }

    #[test]
    fn prepare_gt_ag_forward_no_repeat() {
        let mut f = vec![4u8; 400];
        f[100] = 0;
        f[101] = 0;
        f[102] = 2; // G
        f[103] = 3; // T
        f[198] = 0; // A
        f[199] = 2; // G
        f[200] = 1;
        f[201] = 1;
        let g = make_test_genome(f);
        let pj = prepare_junction(0, 102, 199, 1, &g, g.n_genome);
        assert_eq!(pj.motif, 1); // GT/AG
        assert_eq!(pj.shift_left, 0);
        assert_eq!(pj.shift_right, 0);
        assert_eq!(pj.start_pos, 102);
        assert_eq!(pj.end_pos, 199);
        assert_eq!(pj.strand, 1);
    }

    #[test]
    fn prepare_dot_strand_derived_from_motif() {
        // Non-canonical motif with db_strand=0 → strand 0.
        let g = make_test_genome(vec![0u8; 300]);
        let pj = prepare_junction(0, 100, 200, 0, &g, g.n_genome);
        assert_eq!(pj.motif, 0);
        assert_eq!(pj.strand, 0);

        // GT/AG motif (forward) with db_strand=0 → strand 1.
        let mut f = vec![4u8; 300];
        f[100] = 2;
        f[101] = 3;
        f[199] = 0;
        f[200] = 2;
        let g = make_test_genome(f);
        let pj = prepare_junction(0, 100, 200, 0, &g, g.n_genome);
        assert_eq!(pj.motif, 1);
        assert_eq!(pj.strand, 1);

        // CT/AC motif (reverse) with db_strand=0 → strand 2.
        let mut f = vec![4u8; 300];
        f[100] = 1;
        f[101] = 3;
        f[199] = 0;
        f[200] = 1;
        let g = make_test_genome(f);
        let pj = prepare_junction(0, 100, 200, 0, &g, g.n_genome);
        assert_eq!(pj.motif, 2);
        assert_eq!(pj.strand, 2);
    }

    #[test]
    fn prepare_shift_left_applied_to_coords() {
        // One base of repeat to the left → shift_left=1, coords decrement by 1.
        let mut f = vec![4u8; 400];
        f[102] = 2;
        f[103] = 3;
        f[198] = 0;
        f[199] = 2;
        f[101] = 2; // repeat
        f[100] = 3; // break
        let g = make_test_genome(f);
        let pj = prepare_junction(0, 102, 199, 1, &g, g.n_genome);
        assert_eq!(pj.shift_left, 1);
        assert_eq!(pj.start_pos, 101);
        assert_eq!(pj.end_pos, 198);
    }

    fn pj(
        chr_idx: usize,
        start: u64,
        end: u64,
        motif: u8,
        shift_left: u8,
        strand: u8,
    ) -> PreparedJunction {
        PreparedJunction {
            chr_idx,
            start_pos: start,
            end_pos: end,
            motif,
            shift_left,
            shift_right: 0,
            strand,
        }
    }

    #[test]
    fn stored_and_original_coords_differ_for_canonical_only() {
        let canon = pj(0, 100, 200, 1, 3, 1);
        // Canonical stores ORIGINAL (shifted + shift_left).
        assert_eq!(canon.stored_start(), 103);
        assert_eq!(canon.stored_end(), 203);
        assert_eq!(canon.original_start(), 103);
        assert_eq!(canon.original_end(), 203);

        let noncanon = pj(0, 100, 200, 0, 3, 0);
        // Non-canonical stores SHIFTED.
        assert_eq!(noncanon.stored_start(), 100);
        assert_eq!(noncanon.stored_end(), 200);
        // original_* still recovers pre-shift coords.
        assert_eq!(noncanon.original_start(), 103);
        assert_eq!(noncanon.original_end(), 203);
    }

    #[test]
    fn sort_orders_by_stored_coords() {
        let a = pj(0, 100, 200, 1, 0, 1); // stored 100..200
        let b = pj(0, 50, 150, 1, 0, 1); // stored 50..150
        let c = pj(0, 80, 180, 1, 0, 1); // stored 80..180
        let out = sort_and_dedup(vec![a.clone(), b.clone(), c.clone()]);
        assert_eq!(out, vec![b, c, a]);
    }

    #[test]
    fn dedup_prefers_defined_strand_over_undefined() {
        let defined = pj(0, 100, 200, 1, 0, 1);
        let undefined = pj(0, 100, 200, 0, 0, 0);
        let out = sort_and_dedup(vec![defined.clone(), undefined]);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0], defined);
    }

    #[test]
    fn dedup_collapses_two_non_canonical_to_undefined_strand() {
        // Two non-canonical at same stored coords, opposite strands → one
        // entry with strand = 0.
        let a = pj(0, 100, 200, 0, 0, 1);
        let b = pj(0, 100, 200, 0, 0, 2);
        let out = sort_and_dedup(vec![a, b]);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].strand, 0);
    }

    #[test]
    fn dedup_prefers_canonical_over_non_canonical() {
        // Motif=1 (canonical) beats motif=0 (non) at same stored coords.
        let canon = pj(0, 100, 200, 1, 0, 1);
        let noncan = pj(0, 100, 200, 0, 0, 2);
        let out = sort_and_dedup(vec![noncan, canon.clone()]);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0], canon);
    }

    #[test]
    fn dedup_prefers_strand_matching_motif() {
        // motif=1 (GT/AG +) stored on wrong strand (2) vs a competing
        // motif=2 (CT/AC -) on its correct strand. Same stored coords.
        // STAR keeps the one whose strand matches `2 - motif%2`.
        // For motif=1: expected strand = 2 - 1%2 = 1.
        // For motif=2: expected strand = 2 - 2%2 = 2.
        let old_wrong = pj(0, 100, 200, 1, 0, 2); // motif 1 wants strand 1, has 2
        let new_right = pj(0, 100, 200, 2, 0, 2); // motif 2 wants strand 2, has 2
        let out = sort_and_dedup(vec![old_wrong, new_right.clone()]);
        assert_eq!(out.len(), 1);
        // `old_wrong` is on wrong strand for its motif — STAR replaces.
        assert_eq!(out[0], new_right);
    }
}
