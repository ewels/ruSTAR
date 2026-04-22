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

/// STAR's `sjdbMotif` encoding (see `source/sjdbPrepare.cpp:37-50`).
///
/// `0` = non-canonical / unknown.
/// Values 1-6 correspond to: GT/AG(+), CT/AC(-), GC/AG(+), CT/GC(-),
/// AT/AC(+), GT/AT(-), where each pair is the dinucleotide at the donor /
/// acceptor side of the intron (base-encoded A=0, C=1, G=2, T=3).
pub const MOTIF_NON_CANONICAL: u8 = 0;
pub const MOTIF_GT_AG: u8 = 1;
pub const MOTIF_CT_AC: u8 = 2;
pub const MOTIF_GC_AG: u8 = 3;
pub const MOTIF_CT_GC: u8 = 4;
pub const MOTIF_AT_AC: u8 = 5;
pub const MOTIF_GT_AT: u8 = 6;

/// Classify the splice motif at an intron's donor/acceptor boundary.
///
/// `s` is the 0-based genome-absolute position of the intron's FIRST base
/// (the base immediately after the donor exon's last base). `e` is the
/// 0-based position of the intron's LAST base (immediately before the
/// acceptor exon's first base). Both must be valid indices into `genome`.
///
/// Returns one of the `MOTIF_*` constants.
pub fn classify_motif(genome: &[u8], s: u64, e: u64) -> u8 {
    let si = s as usize;
    let ei = e as usize;
    if ei + 1 > genome.len() || si + 1 >= genome.len() {
        return MOTIF_NON_CANONICAL;
    }
    let b0 = genome[si]; // intron first base
    let b1 = genome[si + 1];
    let b2 = genome[ei - 1];
    let b3 = genome[ei];

    match (b0, b1, b2, b3) {
        (2, 3, 0, 2) => MOTIF_GT_AG,
        (1, 3, 0, 1) => MOTIF_CT_AC,
        (2, 1, 0, 2) => MOTIF_GC_AG,
        (1, 3, 2, 1) => MOTIF_CT_GC,
        (0, 3, 0, 1) => MOTIF_AT_AC,
        (2, 3, 0, 3) => MOTIF_GT_AT,
        _ => MOTIF_NON_CANONICAL,
    }
}

/// Compute STAR's `(sjdbShiftLeft, sjdbShiftRight)` for an intron whose
/// 0-based donor/acceptor positions are `s` and `e`.
///
/// STAR defines these as the number of bases the intron boundary can shift
/// left / right while preserving the donor/acceptor base identity (i.e. the
/// "repeat" length across the junction), capped at 255 per
/// `sjdbPrepare.cpp:52-73`.
///
/// Stops at genome bounds, on any `N` (base code ≥ 4), or when the count
/// reaches 255. Returns `(shift_left, shift_right)` both in `0..=255`.
pub fn compute_shifts(genome: &[u8], s: u64, e: u64, n_genome_real: u64) -> (u8, u8) {
    let si = s as usize;
    let ei = e as usize;
    let ng = n_genome_real as usize;

    // Left shift: keep moving left while G[s-1-jjL] == G[e-jjL] and both
    // bases are A/C/G/T (< 4).
    let mut jj_l: u8 = 0;
    while (jj_l as usize) < si && jj_l < 255 && ei >= jj_l as usize {
        let a_idx = si - 1 - jj_l as usize;
        let b_idx = ei - jj_l as usize;
        if a_idx >= genome.len() || b_idx >= genome.len() {
            break;
        }
        let a = genome[a_idx];
        let b = genome[b_idx];
        if a != b || a >= 4 {
            break;
        }
        jj_l += 1;
    }

    // Right shift: keep moving right while G[s+jjR] == G[e+1+jjR] (< 4).
    let mut jj_r: u8 = 0;
    while (si + jj_r as usize) < ng && jj_r < 255 {
        let a_idx = si + jj_r as usize;
        let b_idx = ei + 1 + jj_r as usize;
        if a_idx >= genome.len() || b_idx >= genome.len() {
            break;
        }
        let a = genome[a_idx];
        let b = genome[b_idx];
        if a != b || a >= 4 {
            break;
        }
        jj_r += 1;
    }

    (jj_l, jj_r)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Build a helper genome buffer. Positions beyond the written slice
    // default to the placeholder padding byte (5) so motif checks at edges
    // are exercised safely.
    fn make_genome(bases: &[(usize, u8)]) -> Vec<u8> {
        let max_idx = bases.iter().map(|(i, _)| *i).max().unwrap_or(0);
        let mut g = vec![5u8; max_idx + 10];
        for (i, b) in bases {
            g[*i] = *b;
        }
        g
    }

    #[test]
    fn classify_gt_ag_forward() {
        // S..E intron = G T ... A G → motif 1
        let g = make_genome(&[(100, 2), (101, 3), (199, 0), (200, 2)]);
        assert_eq!(classify_motif(&g, 100, 200), MOTIF_GT_AG);
    }

    #[test]
    fn classify_ct_ac_reverse() {
        // C T ... A C → motif 2
        let g = make_genome(&[(100, 1), (101, 3), (199, 0), (200, 1)]);
        assert_eq!(classify_motif(&g, 100, 200), MOTIF_CT_AC);
    }

    #[test]
    fn classify_gc_ag() {
        let g = make_genome(&[(100, 2), (101, 1), (199, 0), (200, 2)]);
        assert_eq!(classify_motif(&g, 100, 200), MOTIF_GC_AG);
    }

    #[test]
    fn classify_ct_gc() {
        let g = make_genome(&[(100, 1), (101, 3), (199, 2), (200, 1)]);
        assert_eq!(classify_motif(&g, 100, 200), MOTIF_CT_GC);
    }

    #[test]
    fn classify_at_ac() {
        let g = make_genome(&[(100, 0), (101, 3), (199, 0), (200, 1)]);
        assert_eq!(classify_motif(&g, 100, 200), MOTIF_AT_AC);
    }

    #[test]
    fn classify_gt_at() {
        let g = make_genome(&[(100, 2), (101, 3), (199, 0), (200, 3)]);
        assert_eq!(classify_motif(&g, 100, 200), MOTIF_GT_AT);
    }

    #[test]
    fn classify_non_canonical() {
        // Non-matching combination → 0
        let g = make_genome(&[(100, 0), (101, 0), (199, 0), (200, 0)]);
        assert_eq!(classify_motif(&g, 100, 200), MOTIF_NON_CANONICAL);
    }

    #[test]
    fn classify_out_of_bounds() {
        let g = vec![5u8; 10];
        assert_eq!(classify_motif(&g, 100, 200), MOTIF_NON_CANONICAL);
    }

    fn g_with(bases: &[u8]) -> Vec<u8> {
        // Place `bases` starting at index 100 in a large buffer so the
        // left/right scans have room to run; positions outside get 4 (N)
        // so they break the scan loop naturally.
        let mut g = vec![4u8; 300];
        for (i, b) in bases.iter().enumerate() {
            g[100 + i] = *b;
        }
        g
    }

    #[test]
    fn shifts_no_repeat() {
        // donor side ...GTAG... acceptor — no base match around boundaries.
        // Encoded: A=0 C=1 G=2 T=3. Intron [105,205): GT...AG.
        // Around boundaries: exon1 ends with "A" (0), intron starts with "G" (2), etc.
        let mut g = g_with(&[]);
        g[100] = 0; // donor exon tail
        g[101] = 1;
        g[102] = 0;
        g[103] = 1;
        g[104] = 0;
        // intron
        g[105] = 2; // G
        g[106] = 3; // T
        g[203] = 0; // A
        g[204] = 2; // G
        // acceptor exon head
        g[205] = 3;
        g[206] = 2;
        g[207] = 3;
        g[208] = 2;
        g[209] = 3;
        let (l, r) = compute_shifts(&g, 105, 204, 300);
        assert_eq!(l, 0);
        assert_eq!(r, 0);
    }

    #[test]
    fn shifts_with_repeat_on_left() {
        // Construct a junction where the last base of the donor exon equals
        // the last base of the intron — shifting left by 1 preserves motif.
        // Intron [105, 205). G[104] = X == G[204]; G[103] = Y == G[203].
        let mut g = g_with(&[]);
        // exon1 tail: ... A C
        g[103] = 0; // A
        g[104] = 1; // C
        // intron: G T ... A C (so G[104]=C==G[204], G[103]=A==G[203])
        g[105] = 2;
        g[106] = 3;
        g[203] = 0; // A (matches g[103])
        g[204] = 1; // C (matches g[104])
        // exon2 head
        g[205] = 3;
        let (l, _r) = compute_shifts(&g, 105, 204, 300);
        assert_eq!(l, 2);
    }

    #[test]
    fn shifts_with_repeat_on_right() {
        // STAR check: G[s + jjR] == G[e + 1 + jjR]. So the "first intron
        // base" must equal the "first post-acceptor base". Here s=105,
        // e=204: G[105]==G[205], G[106]==G[206], then the third base
        // differs and the scan stops — shift_right == 2.
        let mut g = g_with(&[]);
        g[105] = 2; // G: intron first base
        g[106] = 3; // T
        g[107] = 0; // A  (will differ from g[207])
        g[204] = 1; // C: intron last base
        g[205] = 2; // G  (matches g[105])
        g[206] = 3; // T  (matches g[106])
        g[207] = 1; // C  (does NOT match g[107])
        let (_l, r) = compute_shifts(&g, 105, 204, 300);
        assert_eq!(r, 2);
    }

    #[test]
    fn shifts_cap_at_255() {
        // Uniform A genome large enough that both scans can run past 255.
        // Scan bounds: left ≤ si; right ≤ genome.len() - (ei + 1).
        // Pick si/ei to leave > 255 bases on each side.
        let g = vec![0u8; 2000];
        let (l, r) = compute_shifts(&g, 500, 1000, 2000);
        assert_eq!(l, 255);
        assert_eq!(r, 255);
    }

    #[test]
    fn shifts_stop_at_n_base() {
        let mut g = vec![0u8; 300];
        // Intron [150, 250). Shift left should stop at N in the scan.
        g[100] = 4; // N: breaks left scan here
        let (l, _) = compute_shifts(&g, 150, 249, 300);
        assert_eq!(l, 49); // scanned 49 A's before hitting N at g[100]
    }
}
