# Experiment 1C: conditional integer-complement diagnosis

## Current status

This work is on branch `codex/full-rank-pppar-datum` in the isolated worktree
`D:\tec\ginan-main-full-rank-ar`.  It starts from commit
`1d91434a9391c529ab87f65b0aae8c4ac814cf50`; the previous Experiment 1
worktree is not modified.

As of 2026-08-31, this is a diagnostic implementation, not a completed PPP-AR
solution.  The new path is disabled by default and never submits second-stage
rows to the Kalman filter.  Runtime Experiment 1C is waiting because another,
unrelated 180-station PEA process is active in WSL; it must not be run
concurrently.

## Independent diagnosis

The receiver datum transform already creates one satellite-minus-pivot integer
coordinate for every non-pivot satellite in each receiver/system/signal group.
For two signals and `m` common satellites, its target dimension is

\[
2(m-1).
\]

The missing rank is introduced later.  `lambda_search()` decorrelates the full
integer vector, multiplies bootstrapped success probabilities from the most
precise end, and retains only the largest tail that still meets the configured
0.9999 threshold.  The corrected 2024 run shows that this tail lies entirely
in the L1C-minus-L2W wide-lane subspace.

Exact trace audits give:

| Run | submitted epochs | submitted rows | exact rank sum | dependent rows | wide-lane rows | single-signal rows | complete-rank epochs | deficit median (range) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| primary | 199/480 | 4,193 | 4,193 | 0 | 4,193 | 0 | 0 | 60 (38–77) |
| independent restart | 178/420 | 3,433 | 3,433 | 0 | 3,433 | 0 | 0 | 64 (41–76) |

The result rules out duplicate submitted rows as the explanation for the rank
deficit.  It also shows that a high partial-fix epoch count is not evidence of
a complete PPP-AR datum.

## Conditional complement construction

Let `a` and `Q` be the full receiver-differenced float integer vector and its
covariance.  LAMBDA supplies a full unimodular transform `Z`, but the existing
partial search submits only its accepted bottom rows.  Partition

\[
Za = \begin{bmatrix}u\\w\end{bmatrix}, \qquad
ZQZ^T = \begin{bmatrix}Q_{uu}&Q_{uw}\\Q_{wu}&Q_{ww}\end{bmatrix},
\]

where `w` is the accepted first-stage integer vector.  For candidate integer
`w_check`, the remaining integer coordinates are diagnosed using

\[
\hat u_{|w}=\hat u+Q_{uw}Q_{ww}^{-1}(w_{check}-\hat w),
\]

\[
Q_{uu|w}=Q_{uu}-Q_{uw}Q_{ww}^{-1}Q_{wu}.
\]

The complement basis is the unused top rows of the same full `Z`.  This choice
is essential: it remains an integer unimodular basis when stacked with the
accepted tail.  An arbitrary real-valued null-space basis would not establish
integer estimability.

The current probe runs LAMBDA on this conditional complement and records only:

- first-stage accepted row count;
- remaining coordinate count;
- second-stage probe row count and validation status;
- combined independent row count and target rank; and
- second-stage bootstrapped success rate and ratio.

Every probe trace is labelled `action=PROBE_ONLY_NOT_SUBMITTED`.

## 2024 development data

Development remains frozen to 2024-07-17 (DOY 199), 30 s sampling, GPS L1C and
L2W, with stations MARS, GRAZ, MATE, DYNG, SOFI and NICO.  G01 remains excluded
because it is absent from the selected WUM orbit, clock and admitted OSB signal
records.  No Zhang `HOU_OSB_LIKE`, GIM or ATT/OBX product is used.

The 23 checksummed inputs are:

1. `igs20_2303.atx`;
2. `igrf14coeffs.txt`;
3. `WUM0MGXFIN_20241990000_01D_01D_ERP.ERP`;
4. `DE436.1950.2050`;
5. `gpt_25.grd`;
6. `OLOAD_GO.BLQ`;
7. `ALOAD_GO.BLQ`;
8. `opoleloadcoefcmcor.txt`;
9. `igs_satellite_metadata_2203_plus.snx`;
10. `sat_yaw_bias_rate.snx`;
11. `IGS0OPSSNX_20241960000_07D_07D_CRD.SNX`;
12. `BRDC00IGS_R_20241990000_01D_MN.rnx`;
13. `WUM0MGXFIN_20241990000_01D_30S_CLK.CLK`;
14. `WUM0MGXFIN_20241990000_01D_01D_OSB.BIA`;
15. WUM final SP3 for DOY 198;
16. WUM final SP3 for DOY 199;
17. WUM final SP3 for DOY 200;
18–23. the six station RINEX observation files.

Exact paths, sizes and SHA-256 values remain in
`experiment_1_2024_input_manifest.json`.  The two new rank-audit JSON files
record the exact per-epoch evidence derived from the corrected formal traces.

## Development gates

The first runtime step is a short window ending just after the first known
accepted event.  It proceeds only when no other PEA process is active.

The short window passes only if:

1. disabling the new option reproduces the corrected first-stage baseline;
2. enabled traces contain only `PROBE_ONLY_NOT_SUBMITTED` complement events;
3. combined independent rows are at least the first-stage count and never
   exceed the target integer rank;
4. conditional covariances remain finite and positive definite;
5. runtime and memory remain suitable for the full 480-epoch development run.

If the whole-network complement search is intractable, the next implementation
will partition by receiver connected component before searching.  It will not
weaken the success-rate or ratio thresholds merely to obtain more fixed rows.

## Work that remains before application

Even a successful probe is insufficient to submit the complement.  Application
requires:

- a joint multi-stage false-fix probability rule, not independent reuse of the
  same nominal threshold without accounting for multiplicity;
- exact stacked-rank and duplicate-row telemetry;
- inverted-sign and double-application phase-OSB controls;
- deliberate wrong-integer injection;
- independent restart/reinitialisation checks on the full integer basis;
- velocity-propagated coordinate acceptance thresholds; and
- a preregistered 2024 confirmation date selected after the algorithm and
  stopping rules are frozen.

No fixed-STEC, full PPP-AR, or scientific covariance claim is made at this
stage.
