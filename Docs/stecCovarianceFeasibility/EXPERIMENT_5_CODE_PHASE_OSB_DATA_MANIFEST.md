# Experiment 5: complete input manifest

All 23 original inputs were hash-verified before the four runs. The audit also checks them after execution. Original paths are under /home/rx/GINAN/inputData in Ubuntu-22.04-G. The two BIA fixtures below replace the original BIA at runtime; all other 22 files remain identical.

| File | Bytes | SHA-256 |
|---|---:|---|
| products/igs20_2303.atx | 40574425 | 5bbc3979e8eae4baa7ea9193a62f979f93e4d433d312f9236a0a054e9c791593 |
| products/tables/igrf14coeffs.txt | 42411 | 8f8d88403028fc4ee92c4f38d97b46e0a87e2cfc496045b43c9e26c1d6b0903c |
| products/tables/DE436.1950.2050 | 9316736 | 9c47f4bef43e0369518f1e0d58207aebb69eeae9401a55b6d98503987fced6c1 |
| products/tables/gpt_25.grd | 578243 | 986fa67d2990befc755b94d6d966baa83b6668d07672d779ea44286e717bee35 |
| products/tables/OLOAD_GO.BLQ | 744986 | 2918b5f155857dd2981822c47ed24e51c30c854cca20412cb6b9ff0d0418fad1 |
| products/tables/ALOAD_GO.BLQ | 238218 | 9a4f2b0a2a3d34ecb3f01e05d2461f54d2a8a65bc65f5d6b6ae42451247a89e9 |
| products/tables/opoleloadcoefcmcor.txt | 24106409 | 3b1d099df46af0c4a7b5d3fd58db609e1b59579fc4329d29f1150b8cf1375e64 |
| products/tables/igs_satellite_metadata_2203_plus.snx | 267785 | b8860c97037befee6002047c32133110d0669bf11686d1b60754729cfd965984 |
| products/tables/sat_yaw_bias_rate.snx | 67094 | 40139be7b046df6e0f8ecf9f3c83e1d9c613c61f786d3119bacf5b06ec1823fd |
| products/IGS0OPSSNX_20241960000_07D_07D_CRD.SNX | 439748 | 3bf384a89d9589c24659b26b1ce667397ad447ad58cbe7956afbefd80f7881f3 |
| products/BRDC00IGS_R_20241990000_01D_MN.rnx | 11069066 | 39f99d229163b02f7ee6f51bbe57067e8da3d9ef6c181ee6845c13a5a1495862 |
| data/MARS00FRA_R_20241990000_01D_30S_MO.rnx | 23779974 | 7e0af5de03234b41332ce9354b5842cdd0fea344d9591b8995ae21039cf306fb |
| data/GRAZ00AUT_R_20241990000_01D_30S_MO.rnx | 31759987 | edf2241d25abf77476f0326807e402d1f32528f85c9e52dfa0358193c5c359d2 |
| data/MATE00ITA_R_20241990000_01D_30S_MO.rnx | 25729678 | 5d706cb8145b7f48236df090fb4b693e9c6b8d7d094707fec87712b595d7d099 |
| data/DYNG00GRC_R_20241990000_01D_30S_MO.rnx | 43911965 | 87ce687d6dca3a373e71ecd58c58a37f334e795666ded8c3c68a5ae6679ef177 |
| data/SOFI00BGR_R_20241990000_01D_30S_MO.rnx | 31709090 | d1a6bfac0d88869a41cf705b90199e278ea2c999d7982447868309da33cf9680 |
| data/NICO00CYP_R_20241990000_01D_30S_MO.rnx | 36989538 | 1568aac7615df081762e61f6eecd64e06c2fd87ec691d6016990448bc56e8232 |
| products/code_2024199/COD0OPSRAP_20241980000_01D_05M_ORB.SP3 | 1456172 | 59254150598c73da61790b0da632228b82f8f835d6e55d06ebcde236ebbec51a |
| products/code_2024199/COD0OPSRAP_20241990000_01D_05M_ORB.SP3 | 1456172 | 38c4d9e79b3622ba4abad4106e062318d03708685065926fc64e3f40078d28b3 |
| products/code_2024199/COD0OPSRAP_20242000000_01D_05M_ORB.SP3 | 1456172 | 8b774390b03a53c60551de47ed79ed656b56c12d2a0f7ae0865a62fc9ef5d072 |
| products/code_2024199/COD0OPSRAP_20241990000_01D_30S_CLK.CLK | 24018732 | bf5b06227be5e9f547f40fcac12aa83390eb15b4872f68d30a629156cd677323 |
| products/code_2024199/COD0OPSRAP_20241990000_01D_01D_OSB.BIA | 940303 | fb99686fff296688920df31e92720d3f66b45c9991aa2ff42f74e732b1542483 |
| products/code_2024199/COD0OPSRAP_20241990000_01D_01D_ERP.ERP | 755 | 1c11b6a64153a61e364811befa52aa70eacf575b0675b2c1cd78492fb6322885 |

## Experimental BIA fixtures

These contain assumed uncertainty, not CODE-published formal sigma. Only 64 GPS L1C/L2W sigma fields change; means and every other byte remain identical. G01 is excluded in processing.

| Assumed sigma | Effective sigma (m) | Fixture SHA-256 |
|---|---:|---|
| 0.003 m | 0.0029999931479602004 | 61dd174b21f91823608634a05dc9184dd1c68840d49d5ac539c0be1d7e05778c |
| 0.01 m | 0.009999997146031201 | 99db024ed3e857756f0ba68fb47a95ca1d8f75d18770b69dfe5466d07b72209f |

Fixture paths: D:/tec/code-exp5-sigma-20260909-v2/assumed_sigma_3mm.BIA and assumed_sigma_10mm.BIA.

The experiment is 2024-07-17, 30 s, 2880 epochs. Stations: MARS, GRAZ, MATE, DYNG, SOFI, NICO. GPS L1C/L2W only. CODE rapid SP3/CLK/ERP/OSB chain; no WUM or CAS products, no ORBEX file. Coordinate reference: the listed IGS weekly CRD SINEX, with no inferred velocity propagation.
