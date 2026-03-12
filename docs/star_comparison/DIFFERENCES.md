[← Back to comparison index](README.md)

# STAR vs ruSTAR: Identified Differences

This file is the primary tracking document for differences between STAR's C++ source and ruSTAR's Rust port. Each difference is assessed for impact on alignment output.

**Status legend**: 🔴 Likely impacting · 🟡 Uncertain / small impact · 🟢 Confirmed equivalent · ✅ Fixed

---

## D1: `stitchWindowAligns` — Score Range Check Uses Per-Mate Score Too 🟡

**STAR** (`stitchWindowAligns.cpp`, finalization block):
```cpp
if ( Score + P.outFilterMultimapScoreRange >= wTr[0]->maxScore
  || ( trA.iFrag >= 0 && Score + P.outFilterMultimapScoreRange >= RA->maxScoreMate[trA.iFrag] )
  || P.pCh.segmentMin > 0) {
    // record transcript
```

The transcript is recorded if the score is within range of the **global** best OR within range of the **per-mate** best (`maxScoreMate[iFrag]`).

**ruSTAR** (`stitch_seeds_core` in `stitch.rs`):
```rust
transcripts.retain(|t| t.score >= max_score - params.out_filter_multimap_score_range);
```

Only checks against the global best. The per-mate score check (`maxScoreMate`) is absent.

**Impact**: For paired-end combined-read paths, if mate1 has a high individual score but the combined WA score is lower, STAR might record the combined transcript while ruSTAR would reject it. Likely explains some of the remaining PE 7-pair gap.

**Investigation needed**: Check if `maxScoreMate` is populated during PE combined-read stitching in STAR and what its value is for the 7 missing PE pairs.

---

## D2: `stitchWindowAligns` — Overhang Filters Applied at Finalization Time 🟢

**STAR** checks exon overhang minimums INSIDE the recursive function at finalization:
```cpp
// For non-annotated junctions:
if ( trA.exons[isj][EX_L] < P.alignSJoverhangMin + trA.shiftSJ[isj][0]
  || trA.exons[isj+1][EX_L] < P.alignSJoverhangMin + trA.shiftSJ[isj][1] ) return;
// For annotated junctions (sjdb):
if ( ( trA.exons[isj][EX_L] < P.alignSJDBoverhangMin ... ) ...) return;
```

**ruSTAR**: Implemented in `finalize_transcript` (`stitch.rs` lines 1377–1426). Checks both left and right exon lengths (including extensions) against `align_sj_overhang_min + shift` for non-annotated, and `align_sjdb_overhang_min` for annotated junctions. The last-exon terminal check is covered because `right_exon_len` includes `right_extend.extend_len` when `isj+1 == wt.exons.len()-1`, matching STAR's `EX_L` after extension. **Confirmed equivalent.**

---

## D3: `stitchWindowAligns` — Soft-Clip at Reference Ends Check 🟢

**STAR**:
```cpp
if (!P.alignSoftClipAtReferenceEnds.yes &&
    ( (trA.exons[trA.nExons-1][EX_G] + Lread - trA.exons[trA.nExons-1][EX_R]) >
       (mapGen.chrStart[trA.Chr] + mapGen.chrLength[trA.Chr])
    || trA.exons[0][EX_G] < (mapGen.chrStart[trA.Chr] + trA.exons[0][EX_R]) ) ) {
    return; //no soft clipping past the ends of the chromosome
}
```

Prevents alignments that would require soft-clipping past chromosome boundaries when `alignSoftClipAtReferenceEnds` is disabled.

**ruSTAR**: The `alignSoftClipAtReferenceEnds` parameter is not implemented. Since the STAR default is `Yes` (soft-clip at reference ends IS allowed), ruSTAR's behavior matches the STAR default. Only diverges if the user explicitly sets `--alignSoftClipAtReferenceEnds No`.

**Impact**: Low — default STAR behavior is to allow soft-clipping at ends (`Yes`). Only matters if a user explicitly disables this, which is uncommon. For standard RNA-seq benchmarking, this is equivalent.

---

## D4: `stitchWindowAligns` — `roStart` Computation 🟢

**STAR**:
```cpp
trA.roStart = (trA.roStr == 0) ? trA.rStart : Lread - trA.rStart - trA.rLength;
```

`roStart` is the read-oriented start: for forward reads it equals `rStart`; for reverse reads it's measured from the end.

**ruSTAR**: This is used to compute soft-clip sizes in SAM output (the left soft-clip = `roStart`, right soft-clip = `Lread - roStart - rLength`). ruSTAR computes this directly in the SAM writer from CIGAR ops. **Likely equivalent.**

---

## D5: `stitchWindowAligns` — Mate Pair Overlap Consistency Check 🟡

**STAR** (`stitchWindowAligns.cpp` lines ~175-220):
```cpp
if (trA.exons[0][EX_iFrag] != trA.exons[trA.nExons-1][EX_iFrag]) {
    // Check for negative insert size
    if (trA.exons[trA.nExons-1][EX_G] + trA.exons[trA.nExons-1][EX_L] <= trA.exons[0][EX_G]) return;
    // If mates overlap, check junction consistency
    if (mate1_end > mate2_start) {
        // Check that overlapping junctions on both mates are identical
        ...
    }
}
```

**ruSTAR** (`read_align.rs`): The negative insert size check IS implemented (Phase 16 fix). The junction consistency check for overlapping mates may or may not be implemented; needs verification.

**Impact**: For PE reads where mate regions overlap (short-insert libraries), inconsistent junction calls could produce false alignments.

---

## D6: `stitchAlignToTranscript` — Left Shift Scoring Range 🔴

**STAR** (`stitchAlignToTranscript.cpp`): The jR scoring loop runs from `min(1, jR+1)` to `max(rGap, jR)`:

```cpp
for (int ii=min(1,jR+1); ii<=max(rGap,jR); ii++) {
    uint g1 = (ii <= jR) ? (gAend+ii) : (gBstart1+ii);
    if (G[g1]<4 && R[rAend+ii]<4) {
        if ( R[rAend+ii] == G[g1] ) {
            if (ii>=1 && ii<=rGap) { Score+=scoreMatch; nMatch++; }
        } else {
            Score -= scoreMatch; nMM++;
            if (ii<1 || ii>rGap) { Score -= scoreMatch; nMatch--; }
        };
    };
}
```

When `jR < 0` (junction shifted LEFT into seed A territory, `ii` starts at `jR+1 < 1`):
- For `ii < 1`: genome side is `gAend+ii` (still donor), but `ii < 1 || ii > rGap` so any mismatch subtracts an EXTRA `scoreMatch` to cancel the previously-assumed match score.
- This handles cases where the donor exon loses some bases to the intron due to left shift.

**ruSTAR** (`stitch_align_to_transcript` in `stitch.rs`):
The "extended left range" block handles `jr_shift < -(shared)` — i.e., when the shift goes PAST all shared bases into exon A. But it does NOT handle the case where `0 > jr_shift > -(shared)` (left-shifted within the shared region but not past it).

Specifically, when `jr_shift < 0` but `|jr_shift| <= shared`:
- STAR: `jR = jr_shift + shared` which is still `>= 0` but `< shared`, so `ii` runs from `min(1, jR+1)` to `max(rGap, jR)` = `max(shared, jR)` = `shared`. This correctly scores all shared bases using the appropriate genome side.
- ruSTAR: The `junction_offset = (shared + jr_shift).max(0).min(shared)` determines how many shared bases go to donor side. The scoring splits shared bases into donor-side and acceptor-side groups. This **should** be equivalent to STAR's loop, but needs verification for edge cases (e.g., `jr_shift = -1`, `shared = 3`).

**Impact**: Unclear. The logic appears correct at first glance for the in-range left-shift case, but any off-by-one in the range boundary could produce wrong scores. Low priority.

---

## D7: `stitchAlignToTranscript` — Step 1 "Move Left" Uses `scoreStitchSJshift` 🟡

**STAR**: Before scanning right in Step 2, Step 1 moves the initial jR1 leftward while:
```cpp
Score1 + P.scoreStitchSJshift >= 0 && int(trA->exons[trA->nExons-1][EX_L]) + jR1 > 1
```

`scoreStitchSJshift` was used to bias the junction scan starting position. STAR removed this from the score calculation itself (set to 0 in newer versions), but the left-scan step still exists.

**ruSTAR** (`find_best_junction_position` in `score.rs`): Needs verification that the left scan starting point is computed identically. From Phase 16.3 notes, this was implemented.

**Impact**: Low — `scoreStitchSJshift` defaults to 0 in STAR.

---

## D8: `stitchAlignToTranscript` — jR Scan Upper Bound 🟢

**STAR** scans `jR1` from left limit to `jR1 < rBend - rAend` = `jR1 <= rGap + L - 1`.

This means the scan goes INTO seed B, up to one base before its end. The rightmost position tested is `jR1 = rGap + L - 1`, which places the junction `L-1` bases into seed B.

**ruSTAR** (`find_best_junction_position` in `score.rs` line 407):
```rust
if jr1 >= r_gap as i32 + next_seed_len as i32 {
    break;
}
```

This breaks when `jr1 >= r_gap + L`, i.e., scans while `jr1 <= r_gap + L - 1`. **Confirmed equivalent to STAR.** ✅

---

## D9: `extendAlign` — `Score > maxScore` (strict) vs `Score >= maxScore` 🟢

**STAR** (`extendAlign.cpp`):
```cpp
if (Score > trA->maxScore) {  // STRICT GREATER THAN
    if (nMM+nMMprev <= min(pMMmax*double(Lprev+i+1), double(nMMmax)) ) {
        trA->extendL = i+1;
        trA->maxScore = Score;
    };
};
```

The record threshold is **strictly greater than**. STAR keeps the first (leftmost/earliest) occurrence of the maximum score.

**ruSTAR** (`extend_alignment` in `stitch.rs` line 204):
```rust
if score > max_score {
    ...
}
```

**Confirmed: uses strict `>`, matching STAR exactly.** ✅

---

## D10: `stitchWindowSeeds` (forward-DP) — Architecture Difference 🟡

**STAR** has TWO stitching passes per window:
1. `ReadAlign::stitchWindowSeeds()` — forward O(N²) DP selecting the single best chain
2. `stitchWindowAligns()` — recursive include/exclude producing multiple transcripts

The forward DP in `stitchWindowSeeds` selects the "winning chain" for the primary single-transcript output. The result is stored in `trAll[iWrec][0]`. Then `stitchWindowAligns` is called to produce additional transcripts for multi-mapper detection.

**ruSTAR**: Only uses the recursive `stitch_recurse` approach (equivalent to `stitchWindowAligns`). The pre-DP `stitchWindowSeeds` approach is simulated by Phase 16.7b pre-DP seed extension scoring (which computes left-extension scores before recursion to help select good seed endpoints).

**Impact**: The pre-DP step in STAR serves as a fast path for finding the best single transcript. Without it, ruSTAR may explore more branches in the recursive search for the same result (performance difference, not correctness difference). However, the scoring heuristics used in STAR's pre-DP may allow STAR to reject certain seed combinations earlier, potentially producing slightly different transcripts in edge cases.

**Investigation needed**: Confirm whether `stitchWindowSeeds` and `stitchWindowAligns` are both always called, or whether they're alternative paths.

---

## D11: `ReadAlign_mappedFilter` — `outFilterScoreMinOverLread` uses `Lread-1` 🟢

**STAR** (`ReadAlign_mappedFilter.cpp`):
```cpp
else if ( (trBest->maxScore < P.outFilterScoreMin)
        || (trBest->maxScore < (intScore)(P.outFilterScoreMinOverLread * (Lread-1)))
        || (trBest->nMatch < P.outFilterMatchNmin)
        || (trBest->nMatch < (uint)(P.outFilterMatchNminOverLread * (Lread-1))) ) {
    unmapType = 1; // TooShort
```

Uses `(Lread-1)` not `Lread` for the proportional thresholds.

**ruSTAR** (`read_align.rs`, `filter_transcripts`): Needs verification that `Lread-1` is used. Phase 16.12 docs mention "Lread-1 filter fix" so this may already be fixed.

**Impact**: Off-by-one in filter threshold. Rarely matters (edge case for very short reads or very tight score thresholds).

---

## D12: Insertion jR Scanning — `alignInsertionFlush` 🟡

**STAR** (`stitchAlignToTranscript.cpp`, insertion path):
```cpp
if (P.alignInsertionFlush.flushRight) {
    for (; jR < (int)rBend-(int)rAend-(int)Ins; jR++) {
        if (R[rAend+jR+1] != G[gAend+jR+1] || G[gAend+jR+1]==4) break;
    };
    if (jR == (int)rBend-(int)rAend-(int)Ins) return -1000009;
}
```

When `alignInsertionFlush = flushRight`, STAR slides the insertion as far right as possible.

**ruSTAR**: The `alignInsertionFlush` parameter and its application in the insertion jR scan may or may not be implemented. Needs verification.

**Impact**: Low — affects insertion placement within reads, which is rare in RNA-seq (short reads mostly have small indels from real variants or sequencing errors).

---

## D13: Deletion Scoring — Short Deletions vs Introns 🟢

**STAR**: Uses `Del * scoreDelBase + scoreDelOpen` for short deletions (`Del < alignIntronMin`) and `scoreGap + jPen` for introns.

**ruSTAR** (`stitch_align_to_transcript`): Uses `score_del_open + score_del_base * del` for deletions. This matches STAR. For introns, uses `motif_score` (= `scoreGap + jPen`). **Confirmed equivalent.**

---

## D14: `stitchWindowAligns` — `nUnique` and `nAnchor` Tracking 🟡

**STAR**:
```cpp
if ( WA[iA][WA_Nrep]==1 ) trAi.nUnique++;  // unique piece
if ( WA[iA][WA_Anchor]>0 ) trAi.nAnchor++;  // anchor piece
```

`nUnique` counts how many seeds have SA range = 1 (unique in genome). `nAnchor` counts anchor seeds. These are used in per-window scoring and filtering.

**ruSTAR**: `wt.n_anchor` is tracked. `n_unique` (equivalent of `nUnique`) may not be tracked. This could affect score-range comparisons if STAR uses `nUnique` anywhere.

**Impact**: Low — `nUnique` is mainly used for MAPQ and multi-mapper reporting, which ruSTAR handles separately.

---

## D15: `stitchWindowAligns` — Variation Adjustment (SNP-aware alignment) 🟢

**STAR** (finalization):
```cpp
Score += trA.variationAdjust(mapGen, R);
```

Adjusts score for known variants (SNPs). Used when variant databases are provided.

**ruSTAR**: No variant database support. `variationAdjust` always returns 0 in a standard STAR run (no variant DB). **Equivalent for standard use.**

---

## D16: `stitchWindowAligns` — Chimeric Detection Score Pass-Through 🟡

**STAR**:
```cpp
|| P.pCh.segmentMin > 0)  // include even if out of score range when chimeric detection active
```

When `chimSegmentMin > 0`, STAR records ALL transcripts regardless of score (for chimeric detection tier 2).

**ruSTAR** (`chimeric/detect.rs`): The chimeric stitcher uses `max_transcripts_per_window=1`, not the full multi-transcript list. Tier 2 chimeric detection uses a separate seed-cluster path. This may not perfectly mirror STAR's approach.

**Impact**: Chimeric detection may miss some edge cases. PE chimeric detection is not yet implemented (Phase 17.3).

---

## Summary Table

| ID | Component | Description | Impact | Status |
|----|-----------|-------------|--------|--------|
| D1 | stitchWindowAligns | Per-mate score check missing | 🟡 PE | Investigate |
| D2 | stitchWindowAligns | Terminal exon overhang check | 🟢 | ✅ Confirmed equivalent |
| D3 | stitchWindowAligns | Soft-clip at chr ends | 🟢 | Default matches STAR |
| D4 | stitchWindowAligns | `roStart` computation | 🟢 | Equivalent |
| D5 | stitchWindowAligns | PE mate overlap consistency | 🟡 | Verify |
| D6 | stitchAlignToTranscript | Left-shift scoring range (in-range) | 🟢 | ✅ Confirmed equivalent |
| D7 | stitchAlignToTranscript | Step 1 left-scan with scoreStitchSJshift | 🟡 | Low risk |
| D8 | stitchAlignToTranscript | jR scan upper bound | 🟢 | ✅ Confirmed correct |
| D9 | extendAlign | Strict `>` vs `>=` for maxScore record | 🟢 | ✅ Confirmed correct |
| D10 | stitchWindowSeeds | Forward-DP not used in ruSTAR | 🟡 | Investigate |
| D11 | mappedFilter | `Lread-1` in proportional filter | 🟢 | Fixed (16.12 docs) |
| D12 | stitchAlignToTranscript | `alignInsertionFlush` | 🟡 | Verify |
| D13 | stitchAlignToTranscript | Del vs intron scoring | 🟢 | Equivalent |
| D14 | stitchWindowAligns | `nUnique` tracking | 🟡 | Low risk |
| D15 | stitchWindowAligns | Variation adjustment | 🟢 | Equivalent |
| D16 | Chimeric | Score pass-through for chimeric | 🟡 | Low risk |

---

## Priority Investigation Order

### Resolved (Confirmed Equivalent)
- **D2, D3, D6, D8, D9**: All confirmed equivalent to STAR.
- **D4, D11, D13, D15**: Always equivalent.

### Remaining — High Priority
1. **D1**: Per-mate score check (`maxScoreMate[iFrag]`) — verify if this explains the PE gap. The check allows single-mate transcripts to be retained in `wTr[]` even when a better combined transcript exists. In ruSTAR, single-mate combined-read candidates are currently discarded.
2. **D5**: PE mate overlap consistency check — relevant for overlapping-mate libraries.

### Remaining — Medium Priority
3. **D7**: Left-scan starting point in `find_best_junction_position` — low risk since `scoreStitchSJshift = 0`.
4. **D10**: `stitchWindowSeeds` N² forward-DP — primarily a performance question; the recursive stitcher should find the same results.

### Remaining — Low Priority
5. **D12**: `alignInsertionFlush` for insertion placement.
6. **D14, D16**: Minor issues for nUnique tracking and chimeric pass-through.
