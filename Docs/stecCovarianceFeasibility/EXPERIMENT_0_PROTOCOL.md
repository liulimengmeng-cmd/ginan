# Experiment 0: fixed-STEC joint-covariance feasibility

## Purpose

Establish whether an unmodified Ginan `main` processing chain can be extended to emit, for every
epoch, the complete posterior covariance of the active `IONO_STEC` state vector without confusing
an ambiguity-resolution attempt with a verified correct integer solution.

This experiment is a prerequisite check. It does not establish that the complete covariance
improves a regional ionospheric-gradient prediction interval.

## Frozen software baseline

- Upstream: `https://github.com/GeoscienceAustralia/ginan.git`
- Base ref: `origin/main`
- Base commit: `7baa32a` (`v4.1.3` release merge)
- Working branch: `codex/main-stec-covariance-feasibility`

The earlier Zhang `HOU_OSB_LIKE` branch is not the software baseline for this experiment. Its
internal products are not IGS Bias-SINEX/OSB products and its prior fixing statistics cannot be
used as acceptance evidence here.

## Source facts established before implementation

1. Ginan's standard `.STEC` writer exports each `IONO_STEC` estimate and its diagonal variance,
   but not cross-covariances.
2. The Kalman state retains the complete covariance matrix `P`; the STEC principal submatrix is
   therefore available at output time.
3. The first-order `IONO_STEC` state is expressed in TECU. The observation coefficient is
   approximately `40.3e16 / f^2` metres per TECU, so the exported covariance unit is TECU squared.
4. Per-epoch ambiguity resolution may operate on a temporary copied state or on the held filter
   state. The covariance exporter must capture the actual state passed to the AR routine.
5. `*_AFTER_AR_ATTEMPT_*` means only that the AR routine ran before capture.
   `*_AFTER_AR_PSEUDOOBS_SUBMITTED_UNVERIFIED` additionally means that the algorithm returned one
   or more resolved combinations and submitted their pseudo-observations to the filter. Neither
   label proves that the integers were correct or that PPP-AR passed an independent validation.
6. RTS sidecars use the corresponding `RTS_SMOOTHED_POSTERIOR_*` label. Forward and smoothed
   covariance blocks must not be pooled without preserving this distinction.

## Output contract

The optional `outputs.ionstec.output_covariance` switch creates a sidecar file with schema marker
`GINAN_STEC_COVARIANCE_V2`.

Each epoch contains:

- one `META` record with state count, upper-triangle count, writer status, numerical asymmetry,
  posterior-stage label, AR invocation flag, eligible ambiguity count, resolved-combination count,
  pseudo-observation submission flag, AR mode, configured validation thresholds, diagnostic status,
  selected decorrelated count, integer-candidate count, actual bootstrap success, best/second
  squared norms and actual solution ratio;
- one `STATE` record per `IONO_STEC` state, including receiver, satellite, state number, filter
  index, estimate in TECU and diagonal variance in TECU squared;
- one `COV` record per upper-triangular matrix entry in TECU squared.

The default limit is 512 STEC states per epoch. If the limit is exceeded, Ginan writes a
`STATE_LIMIT_EXCEEDED` metadata record and no partial matrix. The experiment configuration must
choose a station network small enough to remain below the limit or explicitly raise the limit
after estimating file size.

## Acceptance gates

### G0 — software

- dedicated serializer tests pass;
- full PEA target builds from the frozen base plus this branch;
- generated configuration documentation exposes all three covariance options.

### G1 — file completeness

For every requested epoch:

- writer status is `OK`;
- state indices are contiguous;
- the number of covariance records is exactly `n(n+1)/2`;
- state-diagonal variances equal covariance-matrix diagonal values.

### G2 — numerical covariance validity

Using `scripts/validate_stec_covariance.py`:

- all estimates and entries are finite;
- the captured full matrix was symmetric within the writer tolerance;
- the reconstructed matrix is positive semidefinite within declared numerical tolerance;
- no epoch is accepted from a `STATE_LIMIT_EXCEEDED`, `DIMENSION_MISMATCH`, `NONFINITE_VALUE` or
  `ASYMMETRY_EXCEEDED` writer status.

### G3 — independent fixing evidence

The covariance sidecar's AR fields are control-flow evidence, not proof of a correct fix. A
separate ambiguity-quality table must still establish, per epoch and receiver:

- attempted, accepted and held ambiguity counts;
- ratio/success diagnostics used by the configured Ginan AR method;
- arc age, cycle-slip/reinitialisation status and reference changes;
- whether the exported state was float, partially constrained or fixed.

No epoch may be labelled `fixed STEC` solely because its posterior stage says `AR_ATTEMPT`.

### G4 — continuity

For the Experiment 0 window, at least 70% of scheduled epochs must contain a complete covariance
block and independently verified usable fixed-STEC states. Continuity must be reported by station,
satellite, elevation band and arc, not only as a global average.

### G5 — validation isolation

Validation receivers must be excluded from:

- Ginan network estimation;
- spatial-model fitting;
- covariance calibration;
- date, satellite and event selection.

One held-out receiver cannot provide a two-dimensional gradient truth. Directional-gradient
validation requires a held-out baseline; two-dimensional reconstruction requires at least three
non-collinear held-out receivers.

## Initial run sequence

1. Synthetic serializer and validator tests.
2. Build and configuration-generation check.
3. A 3-minute float smoke run to prove file plumbing only.
4. A 30-minute AR-enabled smoke run with a small GPS L1/L2 network.
5. One quiet-day and one disturbed-day Experiment 0 run after AR quality evidence is available.

Steps 3 and 4 cannot establish the research hypothesis. They only test the processing and output
chain.

The checked-in float smoke overlay is
`Docs/stecCovarianceFeasibility/experiment_0_float_smoke.yaml`. Run it from the repository root and
inject the local data and result roots explicitly:

```bash
LD_LIBRARY_PATH=/home/rx/.local/boost-1.82/lib \
./bin/pea \
  -y Docs/stecCovarianceFeasibility/experiment_0_float_smoke.yaml \
  -a EXP0_DATA_ROOT:/path/to/inputData \
  -a EXP0_OUTPUT_ROOT:/path/to/exp0-float-output
```

The overlay includes `exampleConfigs/ppp_example.yaml`, restricts processing to GPS, disables AR
explicitly, processes six 30-second epochs, and enables both the standard `.STEC` output and the
full-covariance sidecar. It is a plumbing control, not a fixed-STEC experiment.

The companion `experiment_0_ar_smoke.yaml` uses the same three-station data source, processes 60
epochs and invokes main's `LAMBDA_ALT` mode once per epoch without fix-and-hold. It is intended to
expose whether the unmodified main AR path has eligible ambiguities and submits constraints. Even
when its V2 metadata reports submitted pseudo-observations, independent correctness and continuity
gates remain mandatory.

## Stop conditions

Stop and report the prerequisite as failed if any of the following persists:

- no complete STEC covariance block can be produced;
- matrices are non-finite, materially asymmetric or non-PSD;
- the AR output cannot distinguish attempted from accepted fixing;
- usable fixed-STEC continuity remains below 70%;
- validation receivers leak into estimation or calibration.

If the complete covariance later proves negligible relative to the spatial-process covariance,
that is a valid negative scientific result and the main research claim must be revised.
