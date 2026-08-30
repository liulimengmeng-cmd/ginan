# Experiment 0 2024 blind-analysis amendment

## Why this amendment exists

The first FLOAT pivot freeze was created before any validation-only output:

- freeze path: `/home/rx/GINAN/inputData/experiment0_2024/spatial_model_freeze_e0c614d.json`;
- generated: `2026-08-30T10:39:23.277252Z`;
- SHA-256: `7f9f7b8eb68c0139bf9bf33e1d36a5d63f1db5a286023f32134fc6a65e6cda01`;
- declared `heldout_output_read=false`.

The first quiet validation-only PEA run then used only `BALL`, `PARK` and `STR2`.  It ran from
`2026-08-30T10:47:52.17Z` to `2026-08-30T10:49:12.09Z` and exited with status zero.  Before any
held-out residual, coverage, NLL, interval score or estimated `SD_STATE` value was printed or
calculated, a source-code audit identified defects in the frozen evaluator.  The defects were
structural, not data-driven:

1. within-network LOSO calibration omitted the known test--train cross-covariance from the full
   covariance arm;
2. a failed or rank-deficient fit could give the three arms different calibration rows;
3. the freeze did not cryptographically bind the evaluator, its local parser, the complete YAML
   plan or the PEA binary;
4. evaluation inputs lacked pre/post hash stability checks and the report writer allowed
   replacement of an existing report;
5. the full-versus-diagonal coverage comparison did not retain the direction of improvement.

The audit decision was recorded by `2026-08-30T10:53:07Z`, before any storm validation-only PEA
run.  Automated covariance validators had read the quiet matrices only to check serialization,
symmetry, PSD and exact `D x`/`D C D^T` reconstruction.  Their only exposed results were four
`2880/2880 valid` counts.  A sidecar header/META/DATUM check exposed no `SD_STATE` estimate.  This
is structural inspection, but it means the correction is an **analysis-blind amendment**, not a
claim that no program had opened the first quiet held-out files.

## Consequence and forward rule

The original freeze and the first quiet validation-only output are retained as provenance but are
excluded from scientific scoring.  They must not be overwritten or supplied to `evaluate`.

Before further validation processing:

1. repair and test the evaluator using only source structure and the ten-station model outputs;
2. commit the repaired evaluator and protocol;
3. create a new, uniquely named freeze that binds the repaired code, local parser, runtime rules,
   complete configuration closure, manifests/audits, model outputs, PEA binary and two new output
   paths that do not yet exist;
4. run quiet and storm sequentially into those new paths;
5. validate and hash both runs without exposing estimated STEC values;
6. execute the bound evaluator once into a non-existing report path.

Any later method change creates a new amendment, freeze and report name.  It cannot replace an
existing artifact.  The result remains a FLOAT covariance-pipeline Pivot diagnostic because the
registered PPP-AR fixed-availability gate failed; it cannot confirm the original fixed-STEC H1--H5
hypotheses.
