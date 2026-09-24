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
