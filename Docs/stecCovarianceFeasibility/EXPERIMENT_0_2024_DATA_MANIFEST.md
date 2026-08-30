# Experiment 0 2024 Australian data manifest

## Frozen scope

- Official-main base: `origin/main` at `7baa32a`.
- Working branch: `codex/main-stec-covariance-feasibility`.
- Quiet replication: 2024-05-08 (DOY 129), 00:00:00--23:59:30 GPST.
- Storm replication: 2024-05-11 (DOY 132), 00:00:00--23:59:30 GPST.
- Observations: RINEX 3, nominal 30 s, GPS `C1C/L1C/C2W/L2W`.
- Model stations: `ARMC BATH HOB2 MCHL MOBS STR1 SYDN TID1 WGGA YARR`.
- Completely held out: `STR2 BALL PARK`.
- Machine manifest: `/home/rx/GINAN/inputData/experiment0_2024/input_manifest.json`.
- Machine-manifest schema: `GINAN_EXPERIMENT0_2024_INPUT_MANIFEST_V1`.
- Selected-file count: 53.  Every selected file has an uncompressed byte size and SHA-256
  in the machine manifest.
- Frozen machine-manifest SHA-256: `7ba3f1e1661cc2bf6592b4cdb8e10e310ea4bd7eb25120cd33e4a8ecc318f45a`.
- RINEX-QC JSON SHA-256: `52820198dc6926bd07d5dc7998c3b114bb4e6ee863b4c8d485b06e2fcb2a3506`.
- Split-audit JSON SHA-256: `b7a699acb195e8288aa309dd47bb67beb18fbd2da284ddb0b122071414cbaa3c`.
- Product-audit JSON SHA-256: `f5cd9da9a648b65d7df9bcdd70063e1e557148f7f7d5652cfe3d8306eae29c57`.

The held-out RINEX files are downloaded and hashed for a later validation-only run, but are not
listed in either model configuration's `rnx_inputs`.  The split audit reports zero held-out model
inputs for both dates.

## Why these dates

| Index data | Local file | Size/B | SHA-256 | Frozen evidence |
|---|---|---:|---|---|
| GFZ 3-hour Kp | `selection/space_weather/gfz_kp_2024-05-08_2024-05-12.json` | 1,224 | `f6dc7410c552368fd3cb89bf639d7fa6500132f19efd027d9958db9eff199bb5` | DOY129 max Kp 2.0; DOY132 max Kp 9.0 |
| Kyoto provisional hourly Dst | `selection/space_weather/kyoto_dst_provisional_2024-05.txt` | 3,794 | `f2bd464a586e66eb95dd431c5cd864b2c8e180df89eb2d53522fea68ae97f4f7` | DOY129 min +3 nT; DOY132 min -406 nT |

Sources are the [GFZ Kp data service](https://kp.gfz.de/en/data) and the
[Kyoto provisional May 2024 Dst table](https://wdc.kugi.kyoto-u.ac.jp/dst_provisional/202405/index.html).

## Observation files

All observation objects came from the Geoscience Australia GNSS Data Centre API with
`metadataStatus=valid` and `decompress=true`.  Transient signed object URLs are deliberately absent
from the machine manifest.

| Condition | Station | Plain RINEX file | Size/B | SHA-256 |
|---|---|---|---:|---|
| quiet | ARMC | `ARMC00AUS_R_20241290000_01D_30S_MO.rnx` | 59,382,469 | `cd0d231b8c985cc3e7e3125838dbcc8faac198b80fd1ff50d3f7b8a3d316f7b3` |
| quiet | BATH | `BATH00AUS_S_20241290000_01D_30S_MO.rnx` | 29,556,085 | `63449d22d01e61b84f1740c8b6dec05d533278d1512426ae35b6e854977fcce7` |
| quiet | HOB2 | `HOB200AUS_R_20241290000_01D_30S_MO.rnx` | 53,531,419 | `b07c07944b50482c753bab0810dd636c37fa54f0e8bd7614d57416b6fed70e46` |
| quiet | MCHL | `MCHL00AUS_R_20241290000_01D_30S_MO.rnx` | 37,973,450 | `61c581328ff36c5f94af555c4fc37f71e243d1c754592b6c0d0f53fc78e78eb2` |
| quiet | MOBS | `MOBS00AUS_R_20241290000_01D_30S_MO.rnx` | 52,428,965 | `7c22ba7c7ae215c71382ad437ece6f6f6ac4ad2291ad56f106b585d7876b3e06` |
| quiet | STR1 | `STR100AUS_R_20241290000_01D_30S_MO.rnx` | 53,599,133 | `456c06067b82ae69c4f02200c8dee7e0f66ad4ac91a607946bb543dbe2014469` |
| quiet | SYDN | `SYDN00AUS_R_20241290000_01D_30S_MO.rnx` | 41,532,725 | `b950f13066b6c0fdeffeaf55a68320ea115fe4391d40a9184a260fed3b5d27d0` |
| quiet | TID1 | `TID100AUS_R_20241290000_01D_30S_MO.rnx` | 52,070,961 | `1cbcbc01986937fb4fd48396f6d5e115ff135c171b25b3f2bed74bf53b3bd3bc` |
| quiet | WGGA | `WGGA00AUS_S_20241290000_01D_30S_MO.rnx` | 29,930,841 | `cb685480f0f63e9d597e34537185d38a78ba996379216db5b3abfd24f12cf185` |
| quiet | YARR | `YARR00AUS_R_20241290000_01D_30S_MO.rnx` | 56,678,135 | `a67173dc8a60eee382b9c39050983f898acbe6b5b647a8a7b794aba04c483842` |
| quiet | STR2 (held out) | `STR200AUS_R_20241290000_01D_30S_MO.rnx` | 37,229,598 | `6d7e4f92a506cd33720c764ffaf6f06a1fd51bee7c50a90e7c30c85ae15bd16d` |
| quiet | BALL (held out) | `BALL00AUS_S_20241290000_01D_30S_MO.rnx` | 28,353,171 | `372bbb87c5957f9710d269ee0d2ff267bd3c766a34dbde3319a2ec7657ca2abe` |
| quiet | PARK (held out) | `PARK00AUS_R_20241290000_01D_30S_MO.rnx` | 53,390,901 | `756e2acd491d55f1fd80fb84d9d638015aae4768e0c96572e723a3c410db64f0` |
| storm | ARMC | `ARMC00AUS_R_20241320000_01D_30S_MO.rnx` | 58,900,169 | `b9244b42b2fa9ef19cbd5dbbed2a9ad24a7d9fc2c2e96f145dcdc56e3ac015b2` |
| storm | BATH | `BATH00AUS_S_20241320000_01D_30S_MO.rnx` | 29,428,845 | `e50acb71ca679ec9dbf42e3f097cf94e0bad95838ac28d80f75408933c66dbfc` |
| storm | HOB2 | `HOB200AUS_R_20241320000_01D_30S_MO.rnx` | 53,252,659 | `c5f91c0ae8ad85920ed69b06fe603b9899f93b69fffaea67e341784339d8809b` |
| storm | MCHL | `MCHL00AUS_R_20241320000_01D_30S_MO.rnx` | 37,889,592 | `7c9c51b9fc5d1b7e49aa6097381afdab62a6be83b331171ecc5f4dc860d99faf` |
| storm | MOBS | `MOBS00AUS_R_20241320000_01D_30S_MO.rnx` | 52,172,991 | `6191dfd1ad9b580092047403ec3537f0225b76fc5d9ec9cb078874d303a16a89` |
| storm | STR1 | `STR100AUS_R_20241320000_01D_30S_MO.rnx` | 52,656,945 | `bd14e204b68516f91bf43e41765d782f3478dad343b29b82be83c3a142cdc88d` |
| storm | SYDN | `SYDN00AUS_R_20241320000_01D_30S_MO.rnx` | 41,218,151 | `21be0835cbb1d250f9757b09dfcd60bf3ed85a8c07242b09ef0401e0a60f9676` |
| storm | TID1 | `TID100AUS_R_20241320000_01D_30S_MO.rnx` | 51,710,335 | `463ffc25f913387506e7ad6dc8af0d5029a9090194207ef108043aeb80a543ef` |
| storm | WGGA | `WGGA00AUS_S_20241320000_01D_30S_MO.rnx` | 29,784,756 | `fd4a271beb9f0af350e10259784b4550ea2224b049f6bae775e3c960745bfaf9` |
| storm | YARR | `YARR00AUS_R_20241320000_01D_30S_MO.rnx` | 56,392,893 | `ad658f6798a20d60f1cb093e5c4c9427c2fed7cc86832a7ff050d33b609c9894` |
| storm | STR2 (held out) | `STR200AUS_R_20241320000_01D_30S_MO.rnx` | 37,079,108 | `5d0b51bdfd328a0e631ded545e9ab83590641f43aa402309e945d254b58e4870` |
| storm | BALL (held out) | `BALL00AUS_S_20241320000_01D_30S_MO.rnx` | 28,296,887 | `bc9e7457dd1044ee000e2842d4d643bfb3e5e6a4a3d29048ae6a7ac49a6c2449` |
| storm | PARK (held out) | `PARK00AUS_R_20241320000_01D_30S_MO.rnx` | 52,750,845 | `995c0897065cd71dbcc37ccf22a471e895595a1d88a4e97d866f0d0b97e34964` |

RINEX QC is stored in `/home/rx/GINAN/inputData/experiment0_2024/rinex_qc.json`.
Twenty-five files contain 2,880 unique measurement epochs at 30 s with no gap over 30 s.  Storm
WGGA contains 2,873 unique epochs and a maximum 120 s gap; it passes the pre-frozen requirement of
at least 95% coverage and no gap over 300 s.  No duplicate epochs were found.

## Day-specific navigation, precise products and station coordinates

The entire precise chain is WUM final: ORB + CLK + ERP + BIA.  No CODE/WUM mixing is allowed.

| Condition | Role | Relative path | Size/B | SHA-256 |
|---|---|---|---:|---|
| quiet | broadcast navigation | `quiet/products/BRDC00IGS_R_20241290000_01D_MN.rnx` | 11,194,611 | `a0b4b916a9cfde22f61ff88ca602ba19566be3e08b5d8ec8fea44d91eb32b2d9` |
| quiet | previous-day WUM orbit | `quiet/products/WUM0MGXFIN_20241280000_01D_05M_ORB.SP3` | 2,787,204 | `160dea6d0b5720d8ad8cbf4bc4e995a315fb845114f8ca9302ac6d2bd08f01a6` |
| quiet | current-day WUM orbit | `quiet/products/WUM0MGXFIN_20241290000_01D_05M_ORB.SP3` | 2,810,532 | `ed702ded92a1a7dab1ef03a340fa70c8629a4efc26bdd8ebffe9eecf9e296cb0` |
| quiet | next-day WUM orbit | `quiet/products/WUM0MGXFIN_20241300000_01D_05M_ORB.SP3` | 2,810,532 | `933f2d3897ec9c7be05fb7560b05708aa565e762d707dfe870013a4ead05eeba` |
| quiet | WUM 30 s clock | `quiet/products/WUM0MGXFIN_20241290000_01D_30S_CLK.CLK` | 20,750,213 | `9235e57170302b512192d3624561b2f9c9023679c579ad0f41429251b40b729e` |
| quiet | WUM ERP | `quiet/products/WUM0MGXFIN_20241290000_01D_01D_ERP.ERP` | 409 | `6bb6fda2b75d05b73a77d1cf9ab4be22674bb9aaf8a88e3f175f875601dbc4f2` |
| quiet | WUM daily OSB Bias-SINEX | `quiet/products/WUM0MGXFIN_20241290000_01D_01D_OSB.BIA` | 3,143,992 | `55dc975b3a3da2dac1b83937dac24fab5ec6cae3d30b05a36731056665900c77` |
| quiet | APREF daily SINEX | `quiet/products/AUT23133.SNX` | 130,772,295 | `4ecc9ab9ddd2b147615399ceac6ca6aeb1ad596f8ead43615cb42ec8cc51314b` |
| storm | broadcast navigation | `storm/products/BRDC00IGS_R_20241320000_01D_MN.rnx` | 11,685,598 | `0928f7d27ca66efe61d9ad71a9f24b8206b9e36ec994cf0c02421b844209ed52` |
| storm | previous-day WUM orbit | `storm/products/WUM0MGXFIN_20241310000_01D_05M_ORB.SP3` | 2,787,204 | `6cac3fda5f3107be611775058df803d68ebafe5c43b446db8385a7efe4663660` |
| storm | current-day WUM orbit | `storm/products/WUM0MGXFIN_20241320000_01D_05M_ORB.SP3` | 2,833,860 | `c13e61ab46ec51d5c1173da29f4ca72c986d8eabb27db4fa519aab3608de65d7` |
| storm | next-day WUM orbit | `storm/products/WUM0MGXFIN_20241330000_01D_05M_ORB.SP3` | 2,810,532 | `3d69825b11b2c3b6b7d8a56dd6397018341d5e64d05e4828346e8a8224a20800` |
| storm | WUM 30 s clock | `storm/products/WUM0MGXFIN_20241320000_01D_30S_CLK.CLK` | 20,850,253 | `95bdf59583c7903739850cf94ab5c811b7debb9b4e56b77a2920bebf11ce53e2` |
| storm | WUM ERP | `storm/products/WUM0MGXFIN_20241320000_01D_01D_ERP.ERP` | 409 | `da78dac6d8d62a744d81f796fb1b835dce393e1e8227cf5d2ebdeda3239df99b` |
| storm | WUM daily OSB Bias-SINEX | `storm/products/WUM0MGXFIN_20241320000_01D_01D_OSB.BIA` | 3,444,136 | `8829f5e7f667ad4d492d76eb4b70416cfb589fe09fdd6c2bad750d614ea7c4e7` |
| storm | APREF daily SINEX | `storm/products/AUT23136.SNX` | 129,125,877 | `7d873da23d5b88fd9c4424f68a8d4ba3693c9eff1f08bcbdf57b32fb7a42bb93` |

The WUM archive base is
`ftp://igs.gnsswhu.cn/pub/gnss/products/mgex/{GPS_WEEK}/`.  WUM's Bias-SINEX
declares `APC_MODEL IGS20_2303.ATX`.  The eligible product intersection is frozen as:

- quiet: G02--G09 and G11--G32; G01 and G10 excluded;
- storm: G02--G32; G01 excluded.

Eligibility is the intersection of current-day SP3, full-day CLK and BIA `C1C/C2W/L1C/L2W`, not
the BIA satellite list alone.

The independent full-day scan found two real WUM clock gaps:

- quiet G32 lacks 8 records from 23:56:00 through 23:59:30;
- the storm clock file lacks 10 epochs for every satellite from 03:35:00 through 03:39:30.

The local clock hashes match the frozen manifest, so these are product gaps rather than corrupted
downloads.  The WUM final archive has no alternative 1 s, 5 s, 5 min or second 30 s clock file for
either date.  The configuration restricts orbit and clock sources to `PRECISE`; it must not fill the
gaps with broadcast clocks or another analysis centre.  Missing-product epochs remain in the 2,880
epoch denominator and are marked unusable for the fixed-STEC gate.

## Shared static inputs

| Role | Relative path | Size/B | SHA-256 |
|---|---|---:|---|
| WUM-matched antenna calibration | `shared/products/igs20_2303.atx` | 40,574,425 | `5bbc3979e8eae4baa7ea9193a62f979f93e4d433d312f9236a0a054e9c791593` |
| IGS satellite metadata SINEX | `shared/products/tables/igs_satellite_metadata.snx` | 272,883 | `5cae41919b3ec24b464d9ba9855bbf9abd280e0e25c9a671f6833fac0049e6f3` |
| IGRF14 coefficients | `shared/products/tables/igrf14coeffs.txt` | 42,411 | `8f8d88403028fc4ee92c4f38d97b46e0a87e2cfc496045b43c9e26c1d6b0903c` |
| JPL DE436 ephemeris | `shared/products/tables/DE436.1950.2050` | 9,316,736 | `9c47f4bef43e0369518f1e0d58207aebb69eeae9401a55b6d98503987fced6c1` |
| GPT2 grid | `shared/products/tables/gpt_25.grd` | 578,243 | `986fa67d2990befc755b94d6d966baa83b6668d07672d779ea44286e717bee35` |
| ocean tide loading | `shared/products/tables/OLOAD_GO.BLQ` | 744,986 | `2918b5f155857dd2981822c47ed24e51c30c854cca20412cb6b9ff0d0418fad1` |
| atmospheric tide loading | `shared/products/tables/ALOAD_GO.BLQ` | 238,218 | `9a4f2b0a2a3d34ecb3f01e05d2461f54d2a8a65bc65f5d6b6ae42451247a89e9` |
| ocean-pole loading coefficients | `shared/products/tables/opoleloadcoefcmcor.txt` | 24,106,409 | `3b1d099df46af0c4a7b5d3fd58db609e1b59579fc4329d29f1150b8cf1375e64` |
| satellite yaw model | `shared/products/tables/sat_yaw_bias_rate.snx` | 67,094 | `40139be7b046df6e0f8ecf9f3c83e1d9c613c61f786d3119bacf5b06ec1823fd` |

Both daily APREF SINEX files and both loading BLQ files contain all 13 frozen station identifiers.

## Deliberately unselected data

- CODE `COD0OPSRAP` rapid ORB/CLK/ERP/BIA files were downloaded during provider preflight, but they
  are absent from `input_manifest.json` and every configuration.  They are not used in results.
- WUM 30 s BIA is not used; the selected OSB is the daily-constant `01D_01D_OSB.BIA`.
- WUM attitude/ORBEX is not used because the existing ORBEX closure has not been independently
  validated.  Deterministic attitude is retained and eclipse arcs remain a stated limitation.
- No GIM is loaded; first-order slant ionosphere is estimated as a receiver-satellite state.
- No Zhang `HOU_OSB_LIKE` file or Zhang server output is used.  Those are not IGS Bias-SINEX data.

Bias-SINEX v1.00 defines bias as observed minus true, so raw observations are debiased by
subtracting the product bias.  The file has no separate `BIAS_SIGN` field.  See the
[Bias-SINEX v1.00 specification](https://files.igs.org/pub/data/format/sinex_bias_100.pdf).
