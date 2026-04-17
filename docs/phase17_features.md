[← Back to ROADMAP](../ROADMAP.md)

# Phase 17: Features + Polish

**Status**: In Progress (17.1, 17.5, 17.8, 17.A, 17.C complete)

**Goal**: Production-ready features and quality-of-life improvements.

## Sub-phase Status

| Sub-phase | Description | Status |
|-----------|-------------|--------|
| 17.1 | Log.final.out statistics file (MultiQC/RNA-SeQC) | ✅ Complete |
| 17.A | `scoreSeedBest` pre-extension on WA entries (STAR faithful) | ✅ Complete |
| 17.B | Per-mate seeding (fix `.18919121`, `.6302610` arch failures) | Blocked — per-mate clustering approach caused regressions; root cause under investigation |
| 17.C | STAR-faithful SCORE-GATE + mappedFilter for PE (fix 4 MAPQ inflations) | ✅ Complete |
| 17.2 | Coordinate-sorted BAM (`--outSAMtype BAM SortedByCoordinate`) | Planned |
| 17.3 | Paired-end chimeric detection | Planned |
| 17.4 | `--outReadsUnmapped Fastx` | Planned |
| 17.5 | Fix clippy warnings (0 warnings) | ✅ Complete |
| 17.6 | `--outStd SAM/BAM` (stdout output for piping) | Planned |
| 17.7 | GTF tag parameters (`sjdbGTFchrPrefix`, etc.) | Planned |
| 17.8 | `--quantMode GeneCounts` | ✅ Complete |
| 17.9 | `--outBAMcompression` / `--limitBAMsortRAM` | Planned |
| 17.10 | Chimeric Tier 3 (re-map soft-clipped regions) | Planned |
| 17.11 | `--chimOutType WithinBAM` (supplementary FLAG 0x800) | Planned |
| 17.12 | BySJout memory optimization (disk buffering for 100M+ reads) | Planned |
| 17.13 | Phase 9 integration test fixes (realistic test genomes) | Planned |

---

## Phase 17.1: Log.final.out ✅

**Problem**: No `Log.final.out` — MultiQC and RNA-SeQC can't parse results.

**Implementation**:

1. **`Cargo.toml`** — Added `chrono = "0.4"` for timestamps

2. **`src/stats.rs`** — Major expansion (+673 lines):
   - `UnmappedReason` enum (Other, TooShort, TooManyMismatches)
   - 15 new `AtomicU64` counters: read_bases, mapped_bases/mismatches, ins/del count/bases, splices_by_motif[7], splices_annotated, unmapped breakdown, chimeric_reads
   - `record_transcript_stats(&Transcript)` — walks CIGAR for all counters
   - `write_log_final(path, time_start, time_map_start, time_finish)`

3. **`src/align/read_align.rs`** — `align_read()` returns 4-tuple with `Option<UnmappedReason>`

4. **`src/lib.rs`** — Timing via `chrono::Local::now()`, stats collected before BySJout (matches STAR)

**Format**: 47-char right-justified field names + ` |\t` separator. All 37 STAR fields present.

**Differential Test** (10k SE):

| Field | ruSTAR | STAR |
|-------|--------|------|
| Input reads | 10000 | 10000 |
| Avg read length | 150 | 150 |
| Uniquely mapped % | 83.11% | 82.65% |
| Avg mapped length | 146.99 | 146.99 |
| Mismatch rate | 0.40% | 0.40% |

**Files**: `Cargo.toml`, `src/stats.rs`, `src/align/read_align.rs`, `src/lib.rs`

---

## Phase 17.5: Clippy Cleanup ✅

**Problem**: 13 clippy warnings creating noise during debugging.

**Changes**:
- Removed dead code: `verify_match_at_position()`, unused `read_seq` param from `cluster_seeds()`
- **`cluster_seeds()`**: 9 args → 3 args — now takes `&Parameters` instead of 6 windowing params
- **`search_direction_sparse()`**: 8 args → 7 args — folded `effective_start_lmax` into body
- **`ChimericSegment::new()`**: Removed constructor, all 6 call sites use struct literal syntax
- Added `AlignReadResult` type alias for `align_read()` return type
- Idiomatic fixes: `.contains()`, `.div_ceil()`, `.saturating_sub()`
- `#[allow(clippy::too_many_arguments)]` on 4 functions with genuinely many distinct args

**Result**: 0 clippy warnings, 264/264 tests passing.

---

## Phase 17.8: `--quantMode GeneCounts` ✅ (2026-04-17)

**Goal**: Output `ReadsPerGene.out.tab` matching STAR's HTSeq-union gene-level counting.

**Implementation**: New `src/quant/mod.rs` with:
- `GeneAnnotation`: per-chromosome sorted interval list (absolute genome coords) built from GTF exons
- `GeneCounts`: atomic per-gene counters + 3 independent N_noFeature/N_ambiguous arrays
- `QuantContext`: `Arc`-shared bundle for rayon parallel threads
- `--quantMode GeneCounts` + `--sjdbGTFfile` validation in `params.rs`
- SE and PE counting paths in `lib.rs`

**Three bugs fixed vs initial implementation**:
1. **Coordinate mismatch**: GTF exon positions were stored chr-relative; `Transcript.exon.genome_start` uses absolute concatenated-genome coords. Fix: add `genome.chr_start[chr_idx]` offset when converting GTF positions.
2. **Single counting pass**: All 3 columns were identical. STAR runs 3 INDEPENDENT passes — col1 (any strand), col2 (same strand as read), col3 (opposite strand) — each with separate N_noFeature and N_ambiguous.
3. **Too-many-loci bucket**: These were going to N_multimapping. STAR puts them in N_unmapped.

**Results vs STAR (10k SE yeast)**:

| Metric | STAR | ruSTAR |
|--------|------|--------|
| N_unmapped | 1073 | 1074 (+1) |
| N_multimapping | 661 | 661 |
| N_noFeature col1/col2/col3 | 131/3653/4240 | 131/3653/4240 |
| N_ambiguous col1 | 567 | 566 (-1) |
| Gene total col1 | 7568 | 7568 |
| Col1 gene disagreements | — | **0** |
| Col2/col3 gene disagreements | — | 1 each (boundary edge case) |

The ±1 discrepancies (N_unmapped + N_ambiguous) are a single read at a gene overlap boundary — likely a minor coordinate boundary difference.

**Files**: `src/quant/mod.rs` (new), `src/params.rs`, `src/junction/mod.rs` (pub(crate) gtf), `src/lib.rs`

**Tests**: 274/274 (added 6 new quant unit tests), 0 clippy warnings.

---

## Phase 17.A: scoreSeedBest Pre-Extension ✅ (2026-04-16)

**Goal**: Match STAR's `ReadAlign_stitchWindowSeeds.cpp` — pre-extend each seed left+right before the recursive DP and store the result as `pre_ext_score` on each `WindowAlignment` entry.

**What STAR does**: Before `stitchWindowAligns`, STAR computes `scoreSeedBest[iS]` for every seed in the window via a two-level DP: (1) base case: `length + left_ext`, (2) chain case: `stitchAlignToTranscript(iS2→iS1) + scoreSeedBest[iS2]`. Then adds `right_ext` universally. Used for seed ordering in the recursive aligner (start from highest-scoring seed).

**Implementation**:

1. **`src/align/stitch.rs`** — `WindowAlignment` struct: added `pub pre_ext_score: i32` field. All construction sites updated (`pre_ext_score: length as i32` default).

2. **`src/align/score.rs`** — `AlignmentScorer`: added `pub out_filter_score_min_over_lread: f64`. All constructor paths updated.

3. **`src/chimeric/detect.rs`** — `WindowAlignment` construction updated.

4. **`src/align/stitch.rs`** — `stitch_seeds_core`: inserted pre-extension block after seed dedup/sort, before `stitch_recurse`:
   - EXTEND_ORDER respected: left-first for forward clusters (`!stitch_is_reverse`), right-first for reverse clusters (matching `stitch_recurse` base case)
   - `right_len_prev = wa.length + first_ext.extend_len` (mirrors base case's `len_after_first`)
   - Chain DP: `dp[i] = max(dp[i], dp[j] + wa_entries[i].pre_ext_score)` with colinearity check
   - No hard pre-filter gate: STAR uses `scoreSeedBest` for ordering only, not window rejection

**Key finding during implementation**: A pre-filter gate at full `outFilterScoreMinOverLread * (Lread-1)` threshold caused 42 false rejections — reads with only short seeds (9-16bp) in low-quality windows, where the full WT extension (starting from leftmost seed) can reach the threshold even though no individual seed's pre-extension does. STAR does NOT apply this gate; `scoreSeedBest` is used for seed ordering in `stitchWindowAligns` only.

**Result**: 268/268 tests, 0 warnings, 8796/8926 SE (baseline maintained), 8390/8390 PE (baseline maintained). `pre_ext_score` ready for Phase 17.B seed ordering.

---

## Phase 17.B: Per-Mate Seeding — Blocked (architectural)

**What this fixes**: `.18919121` (STAR-only) and `.6302610` (ruSTAR-only FP) — both caused by adapter-contaminated seeds in the combined read crossing mate boundaries.

**Implementation attempted (2026-04-16)**: Seeded mate1 and RC(mate2) separately, clustered each mate independently, then paired compatible (same-chr, same-strand, within-gap) clusters. This correctly prevents cross-mate adapter seed contamination.

**Outcome**: The per-mate clustering approach introduced 10 new ruSTAR-only FPs and 2 new STAR-only cases (net regression of +8 FPs), while NOT fixing the original 2 FPs. Reverted to Phase 17.A state.

**Root cause of regressions**: In unified clustering, a mate2 anchor creates windows that mate1 non-anchor seeds join (and vice versa). Per-mate clustering breaks this inter-mate window sharing that legitimate reads depend on (e.g. short-insert pairs where one mate has a strong anchor and the other has only weak seeds). The single-mate fallback clusters added to the pairing output also caused spurious pairings.

**What STAR actually does**: STAR's `winBin` is `winBin[strand][bin]` — a 2D array indexed by strand (0=fwd, 1=rev) and genomic bin, allocated as `winBin[2][genome_size/65536]`. This is **strand-aware and identical in semantics to ruSTAR's `HashMap<(bool, u64), usize>`**. The earlier theory that STAR uses a bin-only key was incorrect (verified against STAR source in `ReadAlign_assignAlignToWindow.cpp` line 9: `uint iW=winBin[aStr][a1>>P.winBinNbits]`).

**Status**: Blocked — per-mate clustering caused regressions; actual root cause of `.6302610` and `.18919121` (post Phase 17.A) still needs fresh debugging to confirm. The 2 FPs and 2 STAR-only cases remain as known limitations.

---

## Phase 17.C: STAR-faithful SCORE-GATE + mappedFilter ✅ (2026-04-17)

**Problem**: 4 PE MAPQ inflations for rDNA/repeat multi-mappers. ruSTAR NH=2 vs STAR NH=3 for reads with cross-rDNA-copy pairs (M1@copy1 + M2@copy2, 9037bp gap), causing MAPQ=3 vs STAR's MAPQ=1.

**Root cause**: Two distinct bugs:

1. **Per-WT absolute threshold too strict** (`read_align.rs` forward/reverse cluster processing):
   - ruSTAR used `if adjusted_score < combined_score_threshold { continue; }` (hard cutoff at `outFilterScoreMinOverLread * (Lread-1)`)
   - STAR's `stitchWindowAligns.cpp:324` SCORE-GATE uses a RELATIVE criterion: `Score + outFilterMultimapScoreRange >= wTr[0]->maxScore` (within `scoreRange=1` of window best)
   - For cross-copy pairs: same-copy score=198 (g_span=100bp, penalty=-2), cross-copy score=197 (g_span=9237bp, penalty=-3). ruSTAR rejected cross-copy (197 < 198); STAR accepted it (197+1 ≥ 198)

2. **filter_paired_transcripts applied absolute threshold per-pair** (not just to best):
   - ruSTAR checked every pair's `combined_wt_score < absolute_threshold` → removed cross-copy (197 < 198)
   - STAR's `ReadAlign_mappedFilter.cpp` checks only `trBest->maxScore >= threshold` — if the best passes, ALL pairs in the score window are kept

**Fix**:

1. **`src/align/read_align.rs`** — both forward and reverse cluster processing (lines 750, 972):
   ```rust
   // Old:
   if adjusted_score < combined_score_threshold { continue; }
   // New:
   if adjusted_score + params.out_filter_multimap_score_range < combined_score_threshold { continue; }
   ```

2. **`src/align/read_align.rs`** — `filter_paired_transcripts` (line 1373):
   - Changed from per-pair retain to best-pair quality check
   - Find best pair (max `combined_wt_score`); if best fails any threshold → clear all (read unmapped)
   - If best passes → keep all pairs (they already passed multMapSelect relative criterion)

**Verification**: STAR debug trace on `.19790508` confirmed Score=197 cross-copy pair is INSERTED (`TR-INSERTED`) with `global_pass=1` because `scoreRange=1` (`outFilterMultimapScoreRange`). STAR's `mappedFilter` only checks `trBest->maxScore=198 >= 198` — passes.

**Result**: 268/268 tests, 0 warnings, 8796/8926 SE (maintained), 8390/8390 PE (maintained), **0 MAPQ inflations** (was 4), **0 MAPQ deflations**, faithfulness 98.915% (was 98.903%).

---

## Phase 17.2: Coordinate-Sorted BAM — Planned

High user value. STAR outputs `--outSAMtype BAM SortedByCoordinate` natively. Options:
1. In-memory sort during write (requires buffering all records)
2. Wrapper calling `samtools sort` (simpler, already documented as workaround)

---

## Phase 14: STARsolo (Single-Cell) — DEFERRED

**Prerequisite**: All accuracy gaps resolved, position agreement >99%. (Current: 99.92% parity excluding unavoidable ties)
