# Experiment 1 2024 data manifest

## Scope and status

- Experiment: receiver/signal single-difference GPS PPP-AR and STEC covariance export.
- Observation date: 2024-07-17 (2024 day 199), 00:00:00--00:59:30 GPST.
- Sampling: 30 s, 120 epochs.
- Observation signals admitted by the configuration: GPS L1C and L2W only.
- Estimation stations: MARS, GRAZ, MATE, DYNG, SOFI, NICO.
- Reserved held-out stations: ZIM2, WTZR, POTS. They are present locally and covered by the selected SINEX/BLQ files, but are not inputs to this run.
- Configuration: `Docs/stecCovarianceFeasibility/experiment_1_receiver_sd_ar.yaml`.
- Data root used at runtime: `<EXP1_DATA_ROOT>=/home/rx/GINAN/inputData`.
- Source branch: `codex/main-stec-covariance-feasibility`, copied from official `origin/main` at `7baa32a`.
- Current status: the files below are the frozen planned inputs. A post-run section must identify the subset that PEA reports as actually loaded before results are accepted.

The experiment uses only the generic Ginan PPP/PPP-AR path and the external WUM final Bias-SINEX. It does not use `HOU_OSB_LIKE` or any Zhang product/server path.

## Planned input files

Sizes are bytes. SHA-256 values were computed from the uncompressed local files on 2026-08-30.

| Role | Relative path below `/home/rx/GINAN/inputData` | Size | SHA-256 |
|---|---|---:|---|
| Antenna calibration | `products/igs20_2303.atx` | 40574425 | `5bbc3979e8eae4baa7ea9193a62f979f93e4d433d312f9236a0a054e9c791593` |
| Geomagnetic coefficients | `products/tables/igrf14coeffs.txt` | 42411 | `8f8d88403028fc4ee92c4f38d97b46e0a87e2cfc496045b43c9e26c1d6b0903c` |
| Earth rotation | `products/WUM0MGXFIN_20241990000_01D_01D_ERP.ERP` | 409 | `1ee7639f5cd6f5e00d9e8f7bbc2ea757b29b0d116bbbe170d599cd6d34862456` |
| Planetary ephemeris | `products/tables/DE436.1950.2050` | 9316736 | `9c47f4bef43e0369518f1e0d58207aebb69eeae9401a55b6d98503987fced6c1` |
| GPT2 grid | `products/tables/gpt_25.grd` | 578243 | `986fa67d2990befc755b94d6d966baa83b6668d07672d779ea44286e717bee35` |
| Ocean loading | `products/tables/OLOAD_GO.BLQ` | 744986 | `2918b5f155857dd2981822c47ed24e51c30c854cca20412cb6b9ff0d0418fad1` |
| Atmospheric loading | `products/tables/ALOAD_GO.BLQ` | 238218 | `9a4f2b0a2a3d34ecb3f01e05d2461f54d2a8a65bc65f5d6b6ae42451247a89e9` |
| Ocean-pole loading | `products/tables/opoleloadcoefcmcor.txt` | 24106409 | `3b1d099df46af0c4a7b5d3fd58db609e1b59579fc4329d29f1150b8cf1375e64` |
| Satellite metadata | `products/tables/igs_satellite_metadata_2203_plus.snx` | 267785 | `b8860c97037befee6002047c32133110d0669bf11686d1b60754729cfd965984` |
| GPS yaw parameters | `products/tables/sat_yaw_bias_rate.snx` | 67094 | `40139be7b046df6e0f8ecf9f3c83e1d9c613c61f786d3119bacf5b06ec1823fd` |
| Weekly station coordinates | `products/IGS0OPSSNX_20241960000_07D_07D_CRD.SNX` | 439748 | `3bf384a89d9589c24659b26b1ce667397ad447ad58cbe7956afbefd80f7881f3` |
| Broadcast navigation | `products/BRDC00IGS_R_20241990000_01D_MN.rnx` | 11069066 | `39f99d229163b02f7ee6f51bbe57067e8da3d9ef6c181ee6845c13a5a1495862` |
| WUM final clock, 30 s | `products/WUM0MGXFIN_20241990000_01D_30S_CLK.CLK` | 18848226 | `e3dbb996af12276f3496f5085a517d8727bb1526fc67ea31d832ba7806b0b967` |
| WUM final OSB Bias-SINEX | `products/WUM0MGXFIN_20241990000_01D_01D_OSB.BIA` | 6435176 | `faa8da4a6d5813c13c9272b216c8e147a069af71d8db82baefde05ce31f545fa` |
| WUM final orbit, previous day | `products/WUM0MGXFIN_20241980000_01D_05M_ORB.SP3` | 2553924 | `c71561356851fc9ae77111edbf2c86a47147bad5b51ef5334a7cbc0cc15c6186` |
| WUM final orbit, experiment day | `products/WUM0MGXFIN_20241990000_01D_05M_ORB.SP3` | 2553924 | `84eed9d743ab585c98b8fd2e2a36effd575fb197b4a3c2e0072dc33d2a7b1779` |
| WUM final orbit, next day | `products/WUM0MGXFIN_20242000000_01D_05M_ORB.SP3` | 2553924 | `bc922d6f910521b2ae6e5b1b8393de43d1f0f998331c6a3ad06fb69b67de3c9e` |
| MARS observation RINEX | `data/MARS00FRA_R_20241990000_01D_30S_MO.rnx` | 23779974 | `7e0af5de03234b41332ce9354b5842cdd0fea344d9591b8995ae21039cf306fb` |
| GRAZ observation RINEX | `data/GRAZ00AUT_R_20241990000_01D_30S_MO.rnx` | 31759987 | `edf2241d25abf77476f0326807e402d1f32528f85c9e52dfa0358193c5c359d2` |
| MATE observation RINEX | `data/MATE00ITA_R_20241990000_01D_30S_MO.rnx` | 25729678 | `5d706cb8145b7f48236df090fb4b693e9c6b8d7d094707fec87712b595d7d099` |
| DYNG observation RINEX | `data/DYNG00GRC_R_20241990000_01D_30S_MO.rnx` | 43911965 | `87ce687d6dca3a373e71ecd58c58a37f334e795666ded8c3c68a5ae6679ef177` |
| SOFI observation RINEX | `data/SOFI00BGR_R_20241990000_01D_30S_MO.rnx` | 31709090 | `d1a6bfac0d88869a41cf705b90199e278ea2c999d7982447868309da33cf9680` |
| NICO observation RINEX | `data/NICO00CYP_R_20241990000_01D_30S_MO.rnx` | 36989538 | `1568aac7615df081762e61f6eecd64e06c2fd87ec691d6016990448bc56e8232` |

## Deliberately excluded files

- `products/WUM0MGXFIN_20241990000_01D_30S_ATT.OBX`: not loaded because the available ATT file has not passed an independent ORBEX attitude-closure check. Deterministic attitude modelling is used instead.
- `products/IGS0OPSFIN_20241990000_01D_02H_GIM.INX`: not loaded. The first-order slant ionosphere is estimated as a receiver-satellite Kalman state; importing an external GIM would add an unnecessary external ionosphere input to this experiment.
- Zhang `HOU_OSB_LIKE` products and Zhang server outputs: not used. They are not IGS Bias-SINEX products and are outside this official-main experiment.

## Reserved held-out observations

These files are intentionally not listed under `rnx_inputs` for Experiment 1. They will only be hashed and admitted in a separately frozen validation run.

| Station | Local file | Size | Status |
|---|---|---:|---|
| ZIM2 | `data/ZIM200CHE_R_20241990000_01D_30S_MO.rnx` | 17426524 | reserved, not loaded |
| WTZR | `data/WTZR00DEU_R_20241990000_01D_30S_MO.rnx` | 33257084 | reserved, not loaded |
| POTS | `data/POTS00DEU_R_20241990000_01D_30S_MO.rnx` | 60796630 | reserved, not loaded |

## Runtime load audit

Pending. After the clean executable is run, record:

1. executable commit/hash and dirty-state marker;
2. resolved config and all file-open/load messages;
3. any missing antenna, station coordinate, BLQ, OSB signal, orbit, or clock records;
4. the exact PEA-loaded file set and any planned input that was ignored;
5. output file paths and SHA-256 values.

This manifest proves the frozen local byte inputs. It does not by itself prove that ambiguity fixing is correct or that the resulting STEC covariance is statistically calibrated.
