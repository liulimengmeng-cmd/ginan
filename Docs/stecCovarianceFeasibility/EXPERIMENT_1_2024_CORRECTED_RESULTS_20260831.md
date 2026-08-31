# Experiment 1 corrected 2024 PPP-AR results (2026-08-31)

## Result in one sentence

The three identified ambiguity-feedback implementation defects are corrected and
the five-condition 2024 batch has been rerun without another PEA process.  The
corrected primary run submitted checked integer pseudo-observations in 199/480
epochs and the independent restart was consistent on all 6464 shared implied
satellite-pair differences.  All submitted constraints were nevertheless
L1C/L2W wide-lane-type partial constraints: there were zero single-signal rows
and zero complete-rank epochs.  This is therefore runtime evidence for the
software repair, not a complete PPP-AR solution and not fixed-STEC acceptance.

## Frozen software and receipts

- Official upstream: `https://github.com/GeoscienceAustralia/ginan`.
- Official-main base: `origin/main` commit `7baa32a`.
- Working branch: `codex/main-stec-covariance-feasibility`.
- C++ ambiguity-repair commit: `98add0ec3b840a5fb2258800f4c541a1d811872a`.
- Batch source HEAD: `f6c316d2fa54748dbef8343edf9f34b1ef0c478e`.
- Tracked source diff at batch planning: SHA-256
  `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`
  (empty).  The recorded dirty status contains only untracked build products,
  generated version headers and Python caches.
- Frozen PEA SHA-256:
  `6aa8f383a2fa1007364e95222e2e4cf4fb9c5f93b91cd97221cac2c42e098d9e`.
- ELF build ID: `a77e0de05dfe35812f6a366e314cc4ff4493a00f`.
- Embedded PEA version: `untagged-98add0e`.
- Frozen-binary mode: `0555`.
- Authoritative batch root:
  `/home/rx/GINAN/inputData/outputs/exp1_2024_arfix_v2_noconcurrent_20260831_headf6c316d_bin6aa8f383`.
- Batch execution policy: `normal`; active PEA process set was empty at
  preflight and at every per-run resource gate.  No `nice`/`ionice` concurrency
  exception was used.

Artifact hashes:

| Artifact | SHA-256 |
|---|---|
| `batch.plan.json` | `13e4d8a5a2d44c7ed81f0e79df16939f93ca4e5ad2062e3788e375fb675fe673` |
| `batch.receipt.json` | `a4ac26622ba9b24e086e9794954becccceeb195b9536efbe2152adea6c538cbe` |
| `audit/batch.audit.json` | `a76b0ea6eefe0f2f2587e6ba117a075ba1f86746b2f70f35e734fe7d1e3ca0cd` |
| `audit/restart_integer_comparison.json` | `115ce46f0b559f276e80f66f9af24fac51d66d05d6ff497525b2e5a4cd61d340` |
| `audit/coordinate_integrity_primary_vs_float.json` | `2ec7c9690600b45bcfe7db8046aa6b5795ff42a3c00b799d56232ec3a302da8f` |
| `provenance/pea` | `6aa8f383a2fa1007364e95222e2e4cf4fb9c5f93b91cd97221cac2c42e098d9e` |

The batch receipt status is `SUCCESS`.  A separate read-only rehash audit
returned `PASS`, `inputs_match=true`,
`scientific_assessment_matches_receipt=true`, and
`scientific_status=DEFECT_FIX_RUNTIME_EVIDENCE_COMPLETE`.  The audit was then
published with the explicit `--write-audit` option.  The audit command is
read-only by default after commit `c191a7c`.

An earlier batch under
`exp1_2024_arfix_v1_20260831T0205CST_head98add0e_bin6aa8f383` was internally
consistent but ran concurrently with the pre-existing 180-station PEA job and
used `nice -n 5`.  That violated the frozen protocol's no-concurrency rule and
is not the authoritative formal batch.  The v2 batch above removes that
deviation.

## Corrected defects

1. **LAMBDA candidate ordering.**  The previous search could stop after the
   first encountered candidates and could overwrite equal-distance candidates.
   It did not guarantee the global best and second-best squared norms.  The new
   branch-and-bound search retains distinct ties and the globally best requested
   candidates.  A regression changes a false old ratio of 8.148 to the true
   ratio 2.896, which correctly fails a threshold of 3; an equal-distance case
   returns ratio 1.
2. **Feedback design mapping.**  Integer constraints estimated in transformed
   coordinates are now mapped back with `Z*D` before constructing the original
   ambiguity-state pseudo-observations.  Runtime checks compare the constructed
   design rows with that mapping and fail closed on a mismatch.
3. **Pseudo-observation noise.**  Reusing one noise key previously produced the
   rank-one matrix `1e-8 11^T`.  Each resolved row now uses independent direct
   measurement noise, producing `1e-8 I`; runtime checks verify the complete
   matrix before the filter call.

Telemetry labels the final event
`FILTER_CALL_RETURNED_SUBMITTED_UNVERIFIED`.  `filterKalman` returns no
accept/reject value, so the code deliberately does not call the event an
accepted or correct fix.

## Five-condition results

The feedback, mapped-design, diagonal-noise and returned-submission evidence
counts are equal in every AR condition with submitted rows.

| Condition | Epochs | Submitted-unverified epochs | Feedback/design/noise/submission evidence | Maximum candidates | Maximum submitted rows | Complete-rank epochs | Status |
|---|---:|---:|---:|---:|---:|---:|---|
| FLOAT | 480 | 0 | 0/0/0/0 | 0 | 0 | 0 | valid |
| Primary receiver-SD PPP-AR | 480 | 199 | 199/199/199/199 | 2 | 27 | 0 | valid |
| Phase OSB disabled | 480 | 0 | 0/0/0/0 | 2 | 0 | 0 | valid |
| Undifferenced AR negative control | 480 | 17 | 17/17/17/17 | 2 | 10 | 0 | valid stop rule |
| Independent restart at 00:30 | 420 | 178 | 178/178/178/178 | 2 | 26 | 0 | valid |

The zero phase-OSB-disabled result is positive sensitivity evidence for the
presence of phase OSBs, but it does not test an inverted sign or double
application.  The 17 partial undifferenced submissions cannot be promoted: the
negative control has zero complete-rank epochs, and undifferenced ambiguity
states do not close the required integer datum.

## Independent restart comparison

- Primary: 199 submitted-unverified epochs, 4193 constraint rows.
- Restart: 178 submitted-unverified epochs, 3433 constraint rows.
- Overlapping submitted epochs: 140.
- Shared implied satellite-pair comparisons: 6464.
- Integer mismatches: 0.
- Internal graph conflicts: 0.
- Primary/restart wide-lane rows: 4193/3433.
- Primary/restart single-signal rows: 0/0.

The result is `CONSISTENT_SHARED_FIXES` only for the shared wide-lane graph.  It
does not supply an independent certificate that the integers are correct, and
it does not demonstrate narrow-lane or per-signal fixing.

## Coordinate-integrity audit

The coordinate audit selected primary `STATES/AR`, FLOAT `STATES/PPP`, and the
weekly IGS coordinate SINEX.  It parsed 2880 complete station-epoch vectors from
each trace, with no malformed, incomplete or duplicate retained vector.  Exact
primary/FLOAT common coverage is 2880/2880 (six stations times 480 epochs).

| Metric | Count | Maximum | RMS |
|---|---:|---:|---:|
| Primary AR minus FLOAT, 3D | 2880 | 0.011676 m | 0.002745 m |
| FLOAT minus weekly SINEX, 3D | 2880 | 2.612113 m | 0.453407 m |
| Primary minus weekly SINEX, 3D | 2880 | 2.612113 m | 0.453234 m |
| Primary adjacent-epoch jump, 3D | 2874 | 1.501161 m | 0.052457 m |
| Adjacent jump with a receiver submission in the interval | 997 | 0.028204 m | 0.003397 m |
| AR minus FLOAT at a receiver submission | 997 | 0.011676 m | 0.004666 m |

The largest adjacent jump occurs outside a receiver submission interval and is
dominated by early filter convergence, not by a submitted integer constraint.
The audit uses a static weekly SINEX reference without station-velocity
propagation.  No coordinate thresholds were preregistered, so screening status
is `not_evaluated_no_thresholds_configured`.  These numbers do not constitute a
coordinate-integrity scientific pass.

## Exact 2024 input data

- Date: 2024-07-17, GPS day-of-year 199.
- Sampling: 30 s; formal primary/control window 00:00--03:59:30; restart window
  00:30--03:59:30.
- Stations: MARS, GRAZ, MATE, DYNG, SOFI and NICO.
- Signals admitted to ambiguity fixing: GPS L1C and L2W.
- G01 is excluded because it is absent from the selected WUM day-199 orbit,
  clock and both admitted OSB signal records.
- No Zhang `HOU_OSB_LIKE`, GIM or ATT/OBX product is used.

The exact machine-readable list is
`Docs/stecCovarianceFeasibility/experiment_1_2024_input_manifest.json`.  All 23
files matched both size and SHA-256 immediately before the formal run:

| Role | Relative path under `/home/rx/GINAN/inputData` | Bytes | SHA-256 |
|---|---|---:|---|
| antenna calibration | `products/igs20_2303.atx` | 40574425 | `5bbc3979e8eae4baa7ea9193a62f979f93e4d433d312f9236a0a054e9c791593` |
| geomagnetic coefficients | `products/tables/igrf14coeffs.txt` | 42411 | `8f8d88403028fc4ee92c4f38d97b46e0a87e2cfc496045b43c9e26c1d6b0903c` |
| Earth rotation | `products/WUM0MGXFIN_20241990000_01D_01D_ERP.ERP` | 409 | `1ee7639f5cd6f5e00d9e8f7bbc2ea757b29b0d116bbbe170d599cd6d34862456` |
| planetary ephemeris | `products/tables/DE436.1950.2050` | 9316736 | `9c47f4bef43e0369518f1e0d58207aebb69eeae9401a55b6d98503987fced6c1` |
| GPT2 grid | `products/tables/gpt_25.grd` | 578243 | `986fa67d2990befc755b94d6d966baa83b6668d07672d779ea44286e717bee35` |
| ocean loading | `products/tables/OLOAD_GO.BLQ` | 744986 | `2918b5f155857dd2981822c47ed24e51c30c854cca20412cb6b9ff0d0418fad1` |
| atmospheric loading | `products/tables/ALOAD_GO.BLQ` | 238218 | `9a4f2b0a2a3d34ecb3f01e05d2461f54d2a8a65bc65f5d6b6ae42451247a89e9` |
| ocean-pole loading | `products/tables/opoleloadcoefcmcor.txt` | 24106409 | `3b1d099df46af0c4a7b5d3fd58db609e1b59579fc4329d29f1150b8cf1375e64` |
| satellite metadata | `products/tables/igs_satellite_metadata_2203_plus.snx` | 267785 | `b8860c97037befee6002047c32133110d0669bf11686d1b60754729cfd965984` |
| GPS yaw parameters | `products/tables/sat_yaw_bias_rate.snx` | 67094 | `40139be7b046df6e0f8ecf9f3c83e1d9c613c61f786d3119bacf5b06ec1823fd` |
| weekly station coordinates | `products/IGS0OPSSNX_20241960000_07D_07D_CRD.SNX` | 439748 | `3bf384a89d9589c24659b26b1ce667397ad447ad58cbe7956afbefd80f7881f3` |
| broadcast navigation | `products/BRDC00IGS_R_20241990000_01D_MN.rnx` | 11069066 | `39f99d229163b02f7ee6f51bbe57067e8da3d9ef6c181ee6845c13a5a1495862` |
| WUM final 30 s clock | `products/WUM0MGXFIN_20241990000_01D_30S_CLK.CLK` | 18848226 | `e3dbb996af12276f3496f5085a517d8727bb1526fc67ea31d832ba7806b0b967` |
| WUM final Bias-SINEX OSB | `products/WUM0MGXFIN_20241990000_01D_01D_OSB.BIA` | 6435176 | `faa8da4a6d5813c13c9272b216c8e147a069af71d8db82baefde05ce31f545fa` |
| WUM final previous-day orbit | `products/WUM0MGXFIN_20241980000_01D_05M_ORB.SP3` | 2553924 | `c71561356851fc9ae77111edbf2c86a47147bad5b51ef5334a7cbc0cc15c6186` |
| WUM final experiment-day orbit | `products/WUM0MGXFIN_20241990000_01D_05M_ORB.SP3` | 2553924 | `84eed9d743ab585c98b8fd2e2a36effd575fb197b4a3c2e0072dc33d2a7b1779` |
| WUM final next-day orbit | `products/WUM0MGXFIN_20242000000_01D_05M_ORB.SP3` | 2553924 | `bc922d6f910521b2ae6e5b1b8393de43d1f0f998331c6a3ad06fb69b67de3c9e` |
| MARS RINEX observation | `data/MARS00FRA_R_20241990000_01D_30S_MO.rnx` | 23779974 | `7e0af5de03234b41332ce9354b5842cdd0fea344d9591b8995ae21039cf306fb` |
| GRAZ RINEX observation | `data/GRAZ00AUT_R_20241990000_01D_30S_MO.rnx` | 31759987 | `edf2241d25abf77476f0326807e402d1f32528f85c9e52dfa0358193c5c359d2` |
| MATE RINEX observation | `data/MATE00ITA_R_20241990000_01D_30S_MO.rnx` | 25729678 | `5d706cb8145b7f48236df090fb4b693e9c6b8d7d094707fec87712b595d7d099` |
| DYNG RINEX observation | `data/DYNG00GRC_R_20241990000_01D_30S_MO.rnx` | 43911965 | `87ce687d6dca3a373e71ecd58c58a37f334e795666ded8c3c68a5ae6679ef177` |
| SOFI RINEX observation | `data/SOFI00BGR_R_20241990000_01D_30S_MO.rnx` | 31709090 | `d1a6bfac0d88869a41cf705b90199e278ea2c999d7982447868309da33cf9680` |
| NICO RINEX observation | `data/NICO00CYP_R_20241990000_01D_30S_MO.rnx` | 36989538 | `1568aac7615df081762e61f6eecd64e06c2fd87ec691d6016990448bc56e8232` |

## Verification performed

- Full PEA build after the C++ repair: PASS.
- `gnss_ambiguity_diagnostics_tests`: PASS.
- `stec_covariance_tests`: PASS.
- Corrected batch-driver tests: 18/18 PASS on Windows and WSL.
- Coordinate-integrity audit tests: 9/9 PASS under WSL.
- Five PEA conditions and all four covariance validators per condition: PASS.
- Read-only batch rehash audit: PASS.
- `git diff --check` on the changed paths: PASS.

## Scientific boundary and remaining work

Completed means that the identified candidate-search, `Z*D` feedback mapping
and pseudo-observation-noise defects are fixed and exposed by runtime evidence.
It does **not** mean that a full PPP-AR ambiguity solution has been obtained.

The following remain before any fixed-STEC or PPP-AR scientific acceptance:

1. establish a complete integer datum with genuine per-signal/narrow-lane rows,
   rather than only L1C/L2W wide-lane partial constraints;
2. run inverted-sign and deliberate double-application phase-OSB controls;
3. inject a deliberate wrong integer and demonstrate the expected failure;
4. preregister coordinate thresholds and use a velocity-propagated independent
   reference; and
5. repeat the full-rank/arc/slip and fixed-versus-FLOAT scientific assessment
   after the integer-datum problem is closed.
