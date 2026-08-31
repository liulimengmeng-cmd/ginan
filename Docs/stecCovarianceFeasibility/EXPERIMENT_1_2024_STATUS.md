# Experiment 1 2024 status

## Completed on 2026-08-30

- Created the work on `codex/main-stec-covariance-feasibility`, whose official-main base is `7baa32a`.
- Traced the generic PPP ambiguity states to undifferenced cycle-valued `KF::AMBIGUITY` states and confirmed that the previous official-main AR call passed them directly to LAMBDA.
- Added a per-receiver, per-system and per-signal satellite-minus-reference integer transform, with full covariance propagation and explicit reference/singleton diagnostics.
- Extended the STEC covariance sidecar to `GINAN_STEC_COVARIANCE_V3`, including integer-coordinate and receiver-datum accounting.
- Froze a new 2024 day-199 six-station regional data set and its file hashes in `EXPERIMENT_1_2024_DATA_MANIFEST.md`.
- Added float, undifferenced-AR, phase-OSB-disabled and independent-restart controls.

## Verified software checks

- Full `pea` target built successfully.
- `gnss_ambiguity_diagnostics_tests`: PASS, including shared-pivot off-diagonal covariance propagation.
- `stec_covariance_tests`: PASS.
- Python V2/V3 covariance-validator tests: 6/6 PASS under WSL Python.
- Primary and all four control configurations passed PEA dry-run parsing and sanity checks.

The first dry-run exposed that included YAML vector fields are appended. An early overlay therefore retained the 2019 example observations and products. That overlay was rejected and replaced by a standalone configuration. The accepted dry-run summary contains only:

- MARS, GRAZ, MATE, DYNG, SOFI and NICO 2024 day-199 RINEX;
- WUM day-198/199/200 final SP3, day-199 30 s CLK, day-199 ERP and day-199 OSB Bias-SINEX;
- IGS20_2303 ANTEX, IGS weekly coordinate SINEX, broadcast navigation and the listed loading/model tables.

No 2019 observation/product remains in the accepted resolved input list. GIM, ATT/OBX and Zhang internal products are excluded.

## Superseded runtime evidence (audited 2026-08-31)

The statement above that runtime had not started became stale later on
2026-08-30.  A 480-epoch primary run exists at
`/home/rx/GINAN/inputData/outputs/exp1_2024_receiver_sd_primary_480e_503da57`.
Its trace reports `pea_version: untagged-503da57`.  The executable SHA-256 was
not recorded and the executable path was subsequently overwritten, so that
historical binary hash is `NOT_RECORDED/OVERWRITTEN`.

That run predates the 2026-08-31 ambiguity-resolution corrections.  It reported
205/480 forward epochs with submitted integer pseudo-observations, but all were
partial-rank and 199 of the 205 accepted epochs contained only four-term
wide-lane-like constraints.  The old LAMBDA search could overstate the solution
ratio, and all multi-row integer pseudo-observations used a rank-one noise
matrix instead of independent diagonal noise.  The 205/480 value and its
AR-conditioned state/covariance outputs are therefore invalidated.

The correction and rerun requirements were frozen in
`AMBIGUITY_RESOLUTION_CORRECTION_20260831.md`.  The statement that no corrected
runtime result was yet available is superseded by the completed rerun below;
the prohibition on claiming complete PPP-AR acceptance or fixed STEC remains.

## Corrected no-concurrency rerun completed on 2026-08-31

The five-condition corrected batch completed successfully at
`/home/rx/GINAN/inputData/outputs/exp1_2024_arfix_v2_noconcurrent_20260831_headf6c316d_bin6aa8f383`.
No other PEA process was active at preflight or at any per-run resource gate.
The frozen PEA SHA-256 is
`6aa8f383a2fa1007364e95222e2e4cf4fb9c5f93b91cd97221cac2c42e098d9e`.

The primary submitted checked integer pseudo-observations in 199/480 epochs;
the phase-OSB-disabled control submitted none; the undifferenced negative
control submitted partial rows in 17/480 epochs; and the independent restart
submitted rows in 178/420 epochs.  Restart comparison found zero mismatches in
6464 shared implied satellite-pair differences.

Every submitted row is an L1C/L2W wide-lane-type partial constraint.  There are
zero single-signal rows and zero complete-rank epochs in all five conditions.
The software correction is runtime-verified, but a complete PPP-AR solution and
fixed STEC have not been established.  Full results, coordinate audit, exact
input data and artifact hashes are recorded in
`EXPERIMENT_1_2024_CORRECTED_RESULTS_20260831.md`.
