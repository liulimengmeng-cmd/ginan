# Experiment 4: 2024 CODE rapid product manifest

## Experiment window and observations

- Date: 2024-07-17 (DOY 199, GPS week 2323), 30 s, 2880 epochs.
- System/signals: GPS, L1C/L2W.
- Stations: MARS, GRAZ, MATE, DYNG, SOFI, NICO.
- Observation files:
  - `data/MARS00FRA_R_20241990000_01D_30S_MO.rnx` — `7e0af5de03234b41332ce9354b5842cdd0fea344d9591b8995ae21039cf306fb`
  - `data/GRAZ00AUT_R_20241990000_01D_30S_MO.rnx` — `edf2241d25abf77476f0326807e402d1f32528f85c9e52dfa0358193c5c359d2`
  - `data/MATE00ITA_R_20241990000_01D_30S_MO.rnx` — `5d706cb8145b7f48236df090fb4b693e9c6b8d7d094707fec87712b595d7d099`
  - `data/DYNG00GRC_R_20241990000_01D_30S_MO.rnx` — `87ce687d6dca3a373e71ecd58c58a37f334e795666ded8c3c68a5ae6679ef177`
  - `data/SOFI00BGR_R_20241990000_01D_30S_MO.rnx` — `d1a6bfac0d88869a41cf705b90199e278ea2c999d7982447868309da33cf9680`
  - `data/NICO00CYP_R_20241990000_01D_30S_MO.rnx` — `1568aac7615df081762e61f6eecd64e06c2fd87ec691d6016990448bc56e8232`

## Replaced CODE product chain

All paths below are under `products/code_2024199/`.  The files are CODE early
rapid (`COD0OPSRAP`) products archived in the BKG IGS analysis-centre tree.
They are not CODE final products and must not be presented as an equal-latency
CODE-final versus WUM-final comparison.

| Role | File | Bytes | SHA-256 |
|---|---|---:|---|
| Previous-day orbit | `COD0OPSRAP_20241980000_01D_05M_ORB.SP3` | 1,456,172 | `59254150598c73da61790b0da632228b82f8f835d6e55d06ebcde236ebbec51a` |
| Processing-day orbit | `COD0OPSRAP_20241990000_01D_05M_ORB.SP3` | 1,456,172 | `38c4d9e79b3622ba4abad4106e062318d03708685065926fc64e3f40078d28b3` |
| Next-day orbit | `COD0OPSRAP_20242000000_01D_05M_ORB.SP3` | 1,456,172 | `8b774390b03a53c60551de47ed79ed656b56c12d2a0f7ae0865a62fc9ef5d072` |
| Clock | `COD0OPSRAP_20241990000_01D_30S_CLK.CLK` | 24,018,732 | `bf5b06227be5e9f547f40fcac12aa83390eb15b4872f68d30a629156cd677323` |
| Code/phase OSB | `COD0OPSRAP_20241990000_01D_01D_OSB.BIA` | 940,303 | `fb99686fff296688920df31e92720d3f66b45c9991aa2ff42f74e732b1542483` |
| ERP | `COD0OPSRAP_20241990000_01D_01D_ERP.ERP` | 755 | `1c11b6a64153a61e364811befa52aa70eacf575b0675b2c1cd78492fb6322885` |

The audit found full-day GPS coverage for G01–G32 in the day-199 SP3, CLK and
C1C/C2W/L1C/L2W Bias-SINEX records.  G01 nevertheless remains excluded in this
control so that the satellite policy is identical to the WUM experiment.

## Unchanged auxiliary data

- Broadcast navigation: `BRDC00IGS_R_20241990000_01D_MN.rnx`, SHA-256
  `39f99d229163b02f7ee6f51bbe57067e8da3d9ef6c181ee6845c13a5a1495862`.
- Antenna model: `igs20_2303.atx`, SHA-256
  `5bbc3979e8eae4baa7ea9193a62f979f93e4d433d312f9236a0a054e9c791593`.
- Coordinate reference: `IGS0OPSSNX_20241960000_07D_07D_CRD.SNX`, SHA-256
  `3bf384a89d9589c24659b26b1ce667397ad447ad58cbe7956afbefd80f7881f3`.
- IGRF14, DE436, GPT2.5, OLOAD/ALOAD, ocean-pole loading, satellite metadata
  and yaw-bias-rate tables are unchanged from experiment 1.
- No ORBEX attitude file was used; both centres use the same deterministic
  attitude model in this replacement control.

## Product semantics observed

- Bias-SINEX: `BIAS_MODE ABSOLUTE`, `PARAMETER_SPACING 86400`, `TIME_SYSTEM G`,
  `APC_MODEL IGS20`, GPS satellite-clock references `C1W C2W`.
- All 32 GPS L1C and L2W phase-OSB formal standard deviations are zero.  This
  is recorded as “no usable formal product variance”, not as proof of exact
  phase biases.  The Bias-SINEX specification explicitly permits zero standard
  deviation for values taken from an external source.
- C1C formal standard deviation is 0.0039–0.0077 ns; C2W is zero for all 32
  satellites.
- Cross-satellite L1C/L2W phase estimates are strongly related
  (`r=0.998559`, centred slope `1.619559`), but the product does not supply the
  corresponding cross-signal/cross-satellite covariance matrix.
