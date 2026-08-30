# Experiment 0: 2024 Australian STEC-covariance MVP protocol

## Research question

This protocol follows the OSB feasibility conversation's revised question:

> Does the complete, non-diagonal posterior covariance of PPP-AR-derived STEC, together with
> explicit wrong-fix risk, materially change the calibration of regional ionospheric-gradient
> prediction intervals?

The existing 2024 European run is an engineering pretest.  Its 42.7% accepted-epoch rate is below
the registered 70% Experiment 0 prerequisite and is not an MVP pass.

## Software boundary

- Repository: the clean official-main worktree `D:\tec\ginan-main-stec-cov`.
- Branch: `codex/main-stec-covariance-feasibility`.
- Official-main base: `origin/main` at `7baa32a`.
- No Zhang server, Zhang product or `HOU_OSB_LIKE` input is permitted.
- External ambiguity product: WUM final Bias-SINEX OSB.

## Frozen data design

| Item | Quiet replication | Storm replication |
|---|---|---|
| Date | 2024-05-08, DOY129 | 2024-05-11, DOY132 |
| Space weather | max Kp 2.0; min Dst +3 nT | max Kp 9.0; min Dst -406 nT |
| Sampling | GPS L1C/L2W, 30 s, 2,880 scheduled epochs | same |
| Model stations | ARMC BATH HOB2 MCHL MOBS STR1 SYDN TID1 WGGA YARR | identical |
| Completely held out | STR2 BALL PARK | identical |
| Products | WUM final ORB/CLK/ERP/BIA + GA BRDC/APREF | identical series |

The station roles are frozen before processing.  `STR2` is 70.3 m from model station `STR1`,
`BALL` is 103.1 km from `MOBS`, and `PARK` is 130.3 km from `BATH`.  The three held-out stations
also form a non-collinear regional triangle.  Their observations, ambiguities, receiver biases and
STEC states must not enter the model PEA run, hyperparameter fitting, covariance calibration,
method selection or threshold selection.

The exact 53 selected inputs, sizes, hashes and sources are in
`EXPERIMENT_0_2024_DATA_MANIFEST.md` and the machine manifest.  The run must first pass the split
audit with `heldout_input_count=0`.

## Product-gap policy registered before the full run

The WUM clock contains two verified gaps: quiet G32 lacks the final 8 records, while the storm file
lacks 10 epochs for every satellite at 03:35:00--03:39:30.  No same-series replacement exists.

- Only `PRECISE` orbit and clock sources are enabled.
- No broadcast-clock fallback, cross-analysis-centre fill, temporal interpolation or imputation is
  allowed.
- All 2,880 scheduled epochs remain in the denominator.
- An epoch affected by a required product gap is unusable for the fixed-STEC prerequisite even if
  propagated filter states are written.

## Experiment 0 outputs

For every forward and RTS epoch, retain:

1. the original receiver-satellite `IONO_STEC` vector and complete posterior covariance;
2. PPP-AR attempt, datum, success-rate, ratio and pseudo-observation diagnostics;
3. satellite-difference STEC `d_r^{s,q}`;
4. the explicit sparse transform `D`, datum/reference records and globally propagated
   `D C_STEC D^T`, retaining cross-receiver covariance;
5. standard STEC geometry, IPP positions and slant factors.

The original STEC coordinates are filter parameters, not independently certified absolute
geophysical STEC.  The satellite-difference sidecar is the datum-safe quantity used by the spatial
experiment.

## Pre-registered gates

### G0 — data and isolation

- all 26 RINEX files contain the four required GPS observables;
- each file has at least 95% of 2,880 epochs, no duplicate measurement epoch and no gap over 300 s;
- every selected local file matches its manifest SHA-256;
- model/held-out intersection is empty and each model configuration contains exactly the frozen
  10 stations;
- full-day SP3/CLK/BIA signal membership and clock gaps are reported, not inferred from headers.

### G1 — original and satellite-difference covariance

For every claimed usable epoch:

- record counts and dimensions are exact;
- values are finite and the covariance is symmetric within writer/validator tolerance;
- the minimum eigenvalue is no lower than the scale-aware PSD tolerance;
- every difference row has exactly one `+1` target and one `-1` reference term;
- an independent validator reconstructs both `D x` and `D C D^T` from the original sidecar;
- reference changes remain reconstructible from per-epoch `DATUM` and `TRANSFORM` records.

### G2 — candidate fixed-STEC availability

An epoch counts as a candidate usable fixed epoch only when:

- it is not covered by a required-product gap;
- PPP-AR submitted integer pseudo-observations;
- the original and satellite-difference covariance blocks both pass G1;
- the receiver-datum transform has nonzero integer coordinates and no unexplained accounting gap.

At least 2,016 of 2,880 scheduled forward epochs (70%) must qualify on each replication.  This is
an availability gate, not proof that the integers are correct.

### G3 — integer credibility before the word “fixed” is used scientifically

The candidate gate must be followed by independent restart/reinitialization consistency, phase-OSB
sensitivity, float control, arc/pivot continuity, coordinate integrity and a deliberate wrong-fix
or ±1-cycle control.  A high fix rate, LAMBDA acceptance or valid covariance alone is insufficient.

### G4 — held-out construction

Predictions and hyperparameters are frozen before validation-only processing of `STR2/BALL/PARK`.
At least three non-collinear held-out sites or a defensible short-baseline reference must remain.
The validation audit must report zero held-out rows used for fitting, calibration or selection.

## Follow-on no/diagonal/full experiment

Using the same spatial model and split, compare:

- no measurement covariance;
- diagonal-only STEC covariance;
- complete `D C D^T` covariance.

Baselines are IDW, plane/least squares and a spatial covariance model (Matérn or the frozen Zhang
BLUP/Kriging formulation).  Scores include bias, RMSE, NLL, 95% coverage, interval width and
interval score, with block resampling rather than treating 30 s epochs as independent.

Registered scientific thresholds from the conversation:

- full versus zero: at least 5 percentage points less coverage error or at least 10% lower interval
  score, with no more than 25% interval-width increase;
- full versus diagonal: reject a material benefit if coverage differs by less than 3 percentage
  points and NLL by less than 5%;
- calibrated 95% coverage: 93--97% overall and at least 90% in the worst registered subgroup, with
  width inflation no more than 40%;
- a quiet-to-storm coverage loss over 5 percentage points requires storm/elevation/fix-quality
  conditioning rather than a universal calibration claim.

No Go/Pivot/No-Go conclusion is allowed from a 12-epoch smoke test.  The full quiet and storm runs,
independent difference-covariance reconstruction and held-out comparison are mandatory.
