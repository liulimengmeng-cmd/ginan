# PPP-AR ambiguity-resolution correction (2026-08-31)

## Status

The ambiguity-accepted epochs previously reported for Experiment 0 and
Experiment 1 are invalidated.  They were generated before the two software
corrections below and must not be reused, even as partial validation evidence.

The receiver/signal satellite-difference transform itself remains unchanged:

\[
a_D = D a, \qquad Q_D = D Q D^T, \qquad Z_{original} = Z_D D.
\]

Source review and the existing datum-shift/covariance tests found no error in
those three mappings.

## Defect 1: incomplete LAMBDA candidate search

The previous depth-first search stopped as soon as the number of encountered
candidates exceeded `lambda_set_size`.  Encounter order is not squared-norm
order, so the two retained candidates were not guaranteed to be the global
best and second-best integer solutions.  The resulting ratio could be
artificially high and could pass a ratio threshold that the true second-best
candidate failed.

The correction:

- keeps equal-distance candidates as distinct candidates;
- retains the current best configured number of candidates;
- continues the branch-and-bound search after replacing the current worst
  retained candidate; and
- requires at least two candidates for a ratio-based mode.

The focused regression contains both a false-acceptance example (old reported
ratio 8.148, true ratio 2.896, configured threshold 3) and an equal-distance
tie (true ratio 1).  `gnss_ambiguity_diagnostics_tests` passes with the
correction.

Direct evidence of exposure in the old runs is the impossible metadata value
`integerCandidateCount=3` while `lambda_set_size` was the default 2.  It occurs
in 77 of 109 old quiet accepted epochs and 937 of 1738 old storm accepted
epochs.

## Defect 2: rank-one pseudo-observation noise

Every resolved integer row previously used the same `KF::Z_AMB` noise key.
`KFMeas` merged those keys into one scalar noise source, producing

\[
R = 10^{-8} \mathbf{1}\mathbf{1}^T
\]

for multiple constraints, rather than the intended independent noise

\[
R = 10^{-8} I.
\]

All old accepted epochs had at least three submitted combinations and were
affected.  Their AR-conditioned state and STEC covariance outputs are not
valid.

The correction uses direct per-row measurement noise and assigns a distinct
observation number to each row.  Before filtering, the runtime now verifies
that the constructed matrix is exactly the intended independent diagonal
matrix and refuses submission if the check fails.

## Old runtime evidence that is superseded

- Experiment 0 quiet forward: 109/2880 epochs submitted at least one partial
  integer constraint.
- Experiment 0 storm forward: 1738/2880 epochs submitted at least one partial
  integer constraint.
- Experiment 1 primary forward (`untagged-503da57`): 205/480 epochs submitted
  at least one partial integer constraint.
- No accepted epoch in those runs resolved the complete integer-coordinate
  dimension.  Most submitted rows were four-term wide-lane-like constraints,
  not a full per-signal integer datum.

These are defect-exposure counts, not post-correction fixing results.

## Validation state after the corrected rerun

Completed on 2026-08-31:

- full `pea` rebuild; frozen executable SHA-256
  `6aa8f383a2fa1007364e95222e2e4cf4fb9c5f93b91cd97221cac2c42e098d9e`;
- deterministic false-ratio, equal-distance, datum-transform and diagonal-noise
  regression tests;
- corrected primary, FLOAT, phase-OSB-disabled, undifferenced and independent
  restart runs in a serial batch with no other PEA process;
- mapped-design, independent-noise and returned-submission runtime evidence for
  every submitted epoch;
- restart comparison and coordinate-integrity audit; and
- receipt rehash audit with status `PASS`.

The authoritative batch is
`/home/rx/GINAN/inputData/outputs/exp1_2024_arfix_v2_noconcurrent_20260831_headf6c316d_bin6aa8f383`.
Detailed results and hashes are in
`EXPERIMENT_1_2024_CORRECTED_RESULTS_20260831.md`.

The corrected primary run submitted checked pseudo-observations in 199/480
epochs and the restart in 178/420 epochs.  Shared restart comparisons have zero
integer mismatch, and disabling phase OSBs yields zero submitted epoch.
However, all submitted rows are L1C/L2W wide-lane-type partial constraints;
single-signal rows and complete-rank epochs are both zero.  The result is not a
complete PPP-AR solution.

Still pending:

1. a complete per-signal/narrow-lane integer datum;
2. inverted-sign and double-application phase-OSB controls;
3. deliberate wrong-fix injection;
4. preregistered coordinate thresholds with a velocity-propagated reference;
   and
5. full-rank arc/slip and fixed-versus-FLOAT scientific acceptance.

No fixed-STEC claim is made.
