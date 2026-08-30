# Experiment 0 2024 float-covariance pivot protocol

## Trigger and claim boundary

The registered candidate fixed-STEC availability gate failed before any held-out processing:

- quiet forward: 109 of 2,880 epochs submitted integer pseudo-observations;
- storm forward: 1,738 of 2,880 epochs submitted integer pseudo-observations;
- the registered requirement was at least 2,016 epochs on **each** day.

These counts are upper bounds because integer credibility controls have not passed.  The original
fixed-STEC H1--H5 route therefore cannot proceed as a confirmatory experiment.  The pivot below
tests only whether the formal posterior covariance from a float PPP network improves prediction of
independent, datum-safe float PPP STEC contrasts.  It cannot establish correct ambiguity fixing,
wrong-fix risk, an independent geophysical gradient truth or a PPP-AR gradient benefit.

## Frozen separation

- Model observations remain the registered ten stations: `ARMC BATH HOB2 MCHL MOBS STR1 SYDN
  TID1 WGGA YARR`.
- Validation-only observations are exactly `STR2 BALL PARK`.
- The held-out YAMLs include `experiment_0_2024_common.yaml` directly and repeat only the relevant
  day products.  They do not include a model-day YAML because Ginan appends list-valued inputs.
- Ambiguity resolution is explicitly disabled in both model and validation float runs.
- `validate_experiment0_2024_validation_split.py` must pass before a held-out PEA run.
- The model, geometry, covariance arms and quiet model-only calibration are written to a hashed
  freeze JSON before any held-out output is generated or inspected.

## Response and geometry

The primary response is each sidecar `SD_STATE` satellite-minus-reference contrast.  It is not a
point VTEC observation.  Forward `.STEC` supplies only the matching satellite/receiver ECEF
geometry; estimates and covariance always come from the same-epoch float SD sidecar.

The offline geometry definition is fixed as follows:

- mean Earth radius: 6,371.0 km;
- thin-shell height: 506.7 km;
- modified single-layer mapping function with `alpha=0.9782`;
- fixed tangent origin: 33 degrees south, 145 degrees east;
- east/north coordinates are scaled in 1,000 km;
- rows below the already configured 15-degree elevation mask or with missing current geometry are
  rejected without interpolation.

Ginan documents a default mapping height of 506.7 km, but the current `ionmapf` implementation
combines the value with a metre-valued Earth radius without an explicit conversion.  The offline
experiment uses the physically intended 506.7 km definition and must not be described as a
bit-for-bit reproduction of that function.  The first-order `IONO_STEC` state itself is estimated
as slant TEC and does not depend on this offline mapping choice.

For the fixed regional VTEC plane

`h(x) = [1, east_1000km, north_1000km]`,

the design row for a target satellite `s` and the recorded reference `q` is

`A = F_s h(IPP_s) - F_q h(IPP_q)`.

The response is therefore modelled in its actual satellite-difference coordinate.  Pivot changes
are read per epoch; no common pivot is imposed across stations.

## Covariance arms and calibration

All arms use identical epochs, rows and the same three-parameter mean:

- `none`: ordinary least squares, ignoring the model-network SD covariance;
- `diagonal`: GLS with `diag(Q_SD)`;
- `full`: GLS with the complete `Q_SD`.

For scoring, every arm adds the held-out float SD marginal variance.  Cross-covariance between the
separate model and held-out PEA runs is unknown and is not invented.

Quiet model-station leave-one-site-out predictions provide calibration residuals.  Continuous
one-hour blocks receive equal weight in the NLL objective.  Each arm freezes a non-negative
`kappa` by minimizing

`NLL(residual; variance = kappa^2 V_prediction + Q_validation)`.

Held-out observations never enter this step.  Storm data are scored with the quiet calibration and
are not re-calibrated.

## Scores and resampling

Report uncalibrated and calibrated bias, RMSE, Gaussian NLL, marginal 95% coverage, interval width
and interval score, overall and by held-out site.  Comparisons use 1,000 deterministic bootstrap
replicates that resample whole contiguous one-hour blocks (`seed=20240508`), retaining all
station/satellite rows within a sampled block.

The original numerical thresholds are reported as diagnostic analogues only.  Passing them does
not rescue the failed fixed-STEC gate; failing them does not prove that full covariance is useless
for a future credible PPP-AR product.
