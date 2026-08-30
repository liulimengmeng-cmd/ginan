# Experiment 1: 2024 receiver-datum PPP-AR protocol

## Question

Does converting Ginan's per-receiver, per-system and per-signal undifferenced ambiguity states to satellite-minus-reference integer coordinates allow the official-main generic PPP path to submit defensible ambiguity constraints and export the corresponding complete STEC posterior covariance?

This experiment tests an integer-estimability prerequisite. It does not yet establish calibrated regional ionospheric-gradient prediction intervals.

## Frozen software and data

- Branch: `codex/main-stec-covariance-feasibility`.
- Official-main base: `origin/main` commit `7baa32a`.
- Date: 2024-07-17 (day 199), 30 s observations.
- Estimation stations: MARS, GRAZ, MATE, DYNG, SOFI, NICO.
- Signals: GPS L1C and L2W only.
- Satellite screening: G01 excluded because the selected WUM precise orbit, clock and both admitted phase-OSB signals contain no G01 record.
- External products: WUM final SP3/CLK/ERP/OSB, IGS20_2303 ANTEX, IGS weekly coordinate SINEX, broadcast navigation, and the model/loading tables listed with SHA-256 values in `EXPERIMENT_1_2024_DATA_MANIFEST.md`.
- Held-out stations: ZIM2, WTZR, POTS. They do not enter these estimation/control runs.

The external WUM OSB file is Bias-SINEX. No Zhang `HOU_OSB_LIKE` product or Zhang server path is permitted.

## Primary and control runs

| Run | Configuration | Purpose |
|---|---|---|
| Primary | `experiment_1_receiver_sd_ar.yaml` | receiver/signal satellite single differences plus WUM phase OSBs |
| Float | `experiment_1_float_control.yaml` | quantify covariance and coordinate effects without integer constraints |
| Undifferenced AR | `experiment_1_undifferenced_ar_control.yaml` | reproduce the non-integer-estimable official-main AR input as a negative control |
| Phase-OSB disabled | `experiment_1_phase_osb_disabled_control.yaml` | test whether accepted integers depend on satellite phase-OSB application; the BSX file may still be parsed, but its phase correction is disabled |
| Independent restart | `experiment_1_restart_0030.yaml` | reinitialize the filter at 00:30 and compare the overlapping 00:30--03:59:30 window |

All runs use the same observations, precise orbit/clock products, station coordinates, loading models, elevation mask and stochastic settings unless the table explicitly says otherwise.

The primary and three controls run for 480 epochs (00:00--03:59:30).  The
independent restart runs for 420 epochs (00:30--03:59:30).  The former one-hour
and half-hour defaults were too short: the superseded 480-epoch run did not
first reach the pseudo-observation path until 01:39:30.

## Required evidence

For every epoch and receiver/signal group, retain:

1. eligible undifferenced ambiguity count;
2. integer-coordinate count, reference satellite and reference changes;
3. dropped singleton groups;
4. LAMBDA selected count, bootstrapped success, best/second squared norms and ratio;
5. submitted pseudo-observation count;
6. cycle-slip/reinitialisation and arc-age evidence;
7. complete STEC covariance writer/validator status.

The V3 sidecar accounts for the integer-coordinate dimension and receiver-datum groups. The network trace records each reference satellite. A later audit step must convert trace evidence into an epoch/receiver/signal table; the presence of a pseudo-observation alone is not a correct-fix certificate.

## Acceptance and stop rules

- The primary run must never send undifferenced GPS ambiguity states directly to LAMBDA when `receiver_amb_pivot` is true.
- For every V3 epoch, `integer coordinates + datum groups + dropped singleton groups = eligible ambiguities`.
- Every accepted covariance block must be complete, finite, symmetric within writer tolerance and positive semidefinite within validator tolerance.
- The undifferenced and phase-OSB-disabled controls must not be used to promote a fix. If their apparent acceptance is comparable to the primary run, integer correctness remains unresolved.
- The restart overlap must be assessed independently. Agreement only within one continuously propagated filter is not sufficient.
- No `fixed STEC` label is allowed until coordinate integrity, OSB/sign sensitivity, restart consistency, arc continuity and wrong-fix controls all pass.

## Runtime command pattern

Run from the repository root with a unique output root for each condition:

```bash
export LD_LIBRARY_PATH=/home/rx/.local/boost-1.82/lib
./bin/pea \
  -y Docs/stecCovarianceFeasibility/experiment_1_receiver_sd_ar.yaml \
  -a EXP1_DATA_ROOT:/home/rx/GINAN/inputData \
  -a EXP1_OUTPUT_ROOT:/path/to/unique/output
```

Before starting, verify that another memory-intensive PEA network job is not active. Do not run this experiment concurrently with the existing 180-station Zhang experiments in the constrained WSL environment.
