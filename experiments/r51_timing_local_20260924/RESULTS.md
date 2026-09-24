# Local R51 timing: first resumed AR epoch

Date: 2026-09-24. Checkout: `ee3186cb6662796c505dfda720fd367d9cdd520a`.
Instrumented executable SHA-256: `2c921b774d3938a5542aa4d8eeb767687ddcea06d4d84c1c3e8383ec149a67df`.

The run used Ubuntu-22.04-G under WSL (20 GB configured RAM, 10 logical CPUs, 16 GB swap), the frozen six-config R51 180-station input set, and the 2024-07-17 01:59:30 E29 checkpoint. `run_local.sh` used a private mount namespace to redirect the configured output and checkpoint paths to a unique ext4 directory. The original D: output directories were not modified by this test.

Invocation from Windows PowerShell was `wsl -d Ubuntu-22.04-G -u root -- unshare --mount --propagation private bash /mnt/c/Users/rx/Documents/GINAN/ginan-r51-ratio-only-nis-skip-20260923/experiments/r51_timing_local_20260924/run_local.sh`. The script's `run_root` must be changed to a fresh nonexistent directory before re-running it.

Only the first resumed AR epoch, 2024-07-17 02:00:00, was measured. The 02:00:00 `COMMIT.json` reports `ar:true`, and `ZHANG_AR_SUMMARY` reports 2620 candidates and 2080 newly fixed rows. These counts do not by themselves establish a fixed-rate or product-quality claim. After the first product commit the local process was stopped; the one-hour AR run was **not** completed locally.

| Measurement | Time / count |
| --- | ---: |
| AR_TOTAL inclusive | 955.642 s (15m 55.642s) |
| LAYERED_WL_L1_SOLVER inclusive | 767.247 s |
| R51_WL_BLOCK_SEARCH inclusive | 137.1 s |
| R51_L1_BLOCK_SEARCH inclusive | 187.4 s |
| R51_PRODUCT_CLOSURE_FINAL inclusive | 346.804 s |
| R48_OFFICIAL_PRODUCTS inclusive | 200.358 s |
| R48_HISTORY_SEARCH inclusive | 131.368 s |

The inclusive rows are nested and **must not be summed**. The largest accumulated **exclusive** times were:

| Phase | Exclusive time | Calls |
| --- | ---: | ---: |
| KF_ROOT_FACTOR | 143.2 s | 3 |
| R51_PRODUCT_CLOSURE_FINAL | 143.2 s | 1 |
| AR_TOTAL (unscoped remainder) | 109.0 s | 1 |
| R51_WL_PHYSICAL_COMPATIBILITY | 103.1 s | 68 |
| R51_L1_PHYSICAL_COMPATIBILITY | 103.1 s | 57 |
| R48_HISTORY_SEARCH | 82.1 s | 1 |
| R51_PHYSICAL_AFFINE_IMAGE | 57.2 s | 2 |
| R48_PRODUCT_IMAGE | 49.3 s | 2 |
| R48_OFFICIAL_PRODUCTS | 38.5 s | 1 |
| R51_WL_CONDITION | 30.8 s | 68 |
| R48_HISTORY_SUBSET | 29.3 s | 1 |
| R51_L1_CONDITION | 22.5 s | 57 |

The 128 WL ILS calls collectively consumed less than 0.1 s; the ratio search itself was not the first-epoch bottleneck. The complete process used 19m 25s wall time until the controlled stop, including input/setup and post-AR work; `/usr/bin/time` reported mean CPU use of 138% (approximately 1.38 logical CPUs), not ten fully occupied CPUs.

`AR_CALL` began at 7,550,424 KiB RSS and ended at 9,840,068 KiB RSS. Its process high-water mark reached 19,964,772 KiB (about 19.04 GiB), and the end sample reported 299,996 KiB swapped. Therefore the WSL cap is too close to the peak to use later-epoch local timings as clean CPU benchmarks. The local first epoch does not explain the server's later-epoch slowdown without instrumented later-epoch evidence.

Raw evidence remains under `/home/rx/GINAN/local_timing/r51_20260924_instrumented_01/` in Ubuntu-22.04-G: `runner.log`, `outputs/Network-case-202419900.TRACE`, and `outputs/zhang_internal_products.csv.epochs/2024-07-17_02_00_00_ENHANCED_1/COMMIT.json`.

## Writer boundary follow-up

A second isolated first-epoch run used the same inputs/checkpoint and a binary with only additional writer resource boundaries. Its evidence is under `/home/rx/GINAN/local_timing/r51_20260924_writer_probe_01/`. `AR_TOTAL` was 920.863 s; the committed product bundle's `SHA256SUMS` and `COMMIT.json` were byte-identical to the first run. This is diagnostic timing variability, not a proven speedup.

The sampled process high-water mark stayed at 11.53 GiB after solution rows, covariance construction, component gate, and integer-ledger candidate construction. It rose to 19.04 GiB only between `WRITER_AFTER_LEDGER_CANDIDATES` and `WRITER_AFTER_INTEGER_AUTHORITY`, which contains `ProductIntegerLedger::preflight/commit`. At the latter boundary RSS had fallen to 10.18 GiB and 121 MiB was swapped. The transient peak is therefore inside the integer-ledger transaction interval, not the small product CSV/covariance output. These are boundary observations, not allocation-stack proof of the exact expression.

## Ledger transaction follow-up

The feasibility-only membership variant (`r51_20260924_ledger_feasibility_01`) committed a byte-identical first-epoch product bundle (`SHA256SUMS` and `COMMIT.json` matched the writer probe). `AR_TOTAL` was 917.383 s. Its HWM was 19,854,836 KiB versus 19,960,072 KiB in the writer probe: only 105,236 KiB (0.53%) lower. This is not a material memory fix, and one epoch does not establish a speedup.

A subsequent stage-probe run (`r51_20260924_ledger_stage_probe_01`) again produced a byte-identical bundle. Immediately after both the joint HNF and integer-feasibility check, HWM was still 12,093,124 KiB. It remained unchanged after the active-rank check, then reached 17,067,324 KiB after `preflight()` constructed its whole-ledger identity strings and 19,984,120 KiB after `commit()` revalidated them. The measured transaction interval was about 54 seconds. This localizes the extra transient memory to transaction snapshot/commit handling, not to the exact integer decomposition. The writer's identity construction repeatedly serializes the proof closure of every ledger row into large concatenated strings; the stage markers identify the interval, while allocation-level attribution remains unmeasured.

The row-snapshot variant (`r51_20260924_ledger_row_snapshot_01`, executable SHA-256 `986078c9084b40832dfbbaed0ec66f9a2fe6da83f5a800a079385ec4f56cbc16`) stored immutable ledger rows in the preflight receipt and compared their original serialized identity one row at a time during commit. The 02:00:00 committed `SHA256SUMS` and `COMMIT.json` matched the writer-probe baseline byte-for-byte. `AR_TOTAL` was 891.533 s. HWM remained 12,091,868 KiB through `LEDGER_AFTER_SNAPSHOT_IDENTITY`, `WRITER_AFTER_INTEGER_AUTHORITY`, and writer exit, with zero swap. Relative to the writer-probe HWM of 19,960,072 KiB, this is 7,868,204 KiB (7.50 GiB, 39.4%) less peak resident memory. The 29.3 s shorter `AR_TOTAL` is a single-run observation, not a statistically established speedup. This local run was stopped after the first product commit; it does not establish sustained one-hour memory behavior or product quality.

Validation: the targeted R46 integer-feasibility and stale-preflight tests passed (3 cases, 14 assertions), and the three R51 performance/entailment/compatibility CTests passed. The complete `zhang_full_rank_tests` run passed 450 of 451 cases (17,356 of 17,357 assertions). The sole failure was the previously observed `e29_checkpoint_rejects_corruption_and_provenance_drift` assertion at `test_ZhangFullRank.cpp:9576`: it expects bare `CHECKPOINT_PROVENANCE_MISMATCH`, while the result includes `:CONFIG_SHA256:expected=...:actual=...`. The suite is therefore **not fully green**. No one-hour AR acceptance test was run for this change.
