# R51 physical-arc transport regression

The graph lifecycle repair preserves surviving physical phase functionals when
the old tree cannot be reused. It supports a rectangular state transform and
independent finite priors for fresh arcs, retaining full covariance with other
states. It does not certify new integers or change acceptance/noise thresholds.

The fallback now transports component coordinates before admitting carrier
rows. A disconnected FLOAT forest keeps uncertain component offsets; the
existing single-component AR API rejects it. Only explicit legacy/shadow
controls retain the old node-level reset. Nonconstant phase process models are
rejected until their dynamics can also be transported.

Two performance defects are repaired: a transform no longer resurrects pending
state removals, and terminal reconciliation runs only when the posterior was
solved against an older model generation. Last-iteration model changes still
require reconciliation.

## Reproduce after a Release Unix Makefiles build

Build and run the `zhang_full_rank_tests` target using the repository's existing
test setup. Then, from the repository root, run:

```sh
python3 src/tests/run_r51_algebra_transport.py build
```

The second test links the actual PEA objects with a separate test driver and a
renamed copy of the production main object. It does not replace the production
executable or use a stub KF implementation. It exercises:

- scheduled STEC removal and marginal covariance preservation;
- correlated fresh priors, factor routing, commit sequence and rejection rollback;
- retired-chord audits and rejection of retired replacement tree edges;
- skipping a redundant converged-model solve while retaining final-change reconciliation;
- production graph-controller recovery when all old physical arcs retire together.

The mathematical suite additionally checks simultaneous tree/chord retirement,
new-node augmentation, split and rejoined components, physical-row identities,
cross-covariance preservation, and unsupported dynamic-model rejection.

## Evidence limits

The 2026-09-16 local validation passed 439 test cases and 17,214 assertions.
A 180-station, 20-epoch comparison against the frozen parent program exercised
the component fallback at 2024-07-17 00:08:30 GPST. The repaired run admitted
1397/1397 eligible phase edges and retained 2780 single-frequency physical
relations; the parent used `local_reinitialise` at that epoch. Those counts are
not fixed ambiguities. Both runs had one nonconverged QC epoch.

With BLAS=1 and OMP=4, total wall times were 1153.66 s (repaired) and 1127.56 s
(parent). Thus this comparison does not demonstrate an overall speedup, despite
terminal reconciliations falling from seven to one. Final boundary fixes then
passed the production controller test and a fresh single-epoch smoke test with
the same first-epoch posterior digest.

The final executable subsequently completed all 20 epochs with BLAS=4 and
OMP=4 in 517.45 s, compared with 1127.56 s for the frozen parent with BLAS=1
(54.1% less wall time). Its last 18 mean commit intervals were 16.52 s versus
34.78 s for the parent. It performed five component reparameterizations, no
local resets, one terminal reconciliation and one nonconverged QC epoch.
This comparison changes both code and the BLAS setting; it is not a pure
thread-scaling benchmark or a day-long runtime prediction. The final program
and the repaired BLAS=1 run have the same 1160 product keys and validity flags;
correction differences are at most 0.867 mm (0.110 mm RMS), not truth errors.
All comparison runs disable checkpoint writing. The original 48-hour run
remains stopped.

The full local evidence is under `r51_graph_speed_fix_20260916`, outside this
checkout, including the source conversation, configuration, commands, hashes,
logs, per-epoch comparison, and any subsequent threading experiment. Long-term
state retirement, 24/48-hour stability, physical FLOAT accuracy and correct AR
acceptance require separate validation. Old checkpoints are not relabeled as
compatible with a changed executable.
