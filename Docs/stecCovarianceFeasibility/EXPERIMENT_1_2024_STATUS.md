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

## Pending runtime evidence

The actual 2024 PEA smoke run has not started. At the time of the readiness check, PID 51 was an existing 180-station Zhang PEA run consuming about 8.8--10.2 GB RAM and more than eight CPU cores. Starting another PEA run would risk both experiments, so it was deliberately deferred without terminating or altering the existing process.

After resources are free, the next required actions are:

1. commit and rebuild a clean executable so the run hash has no dirty marker;
2. run a 12-epoch primary smoke and audit the files actually loaded;
3. validate V3 covariance blocks and receiver-datum accounting;
4. run the full primary/control/restart matrix;
5. assemble arc, slip, reference-change and integer-quality evidence before any `fixed STEC` claim.

This status is software/configuration evidence only. It is not PPP-AR acceptance and not evidence that complete covariance improves regional ionospheric-gradient prediction intervals.
