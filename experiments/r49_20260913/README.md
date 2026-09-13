# R49 server experiment

Read the local source-conversation export and R49_IMPLEMENTATION.md for scope.
Configuration has the same 180 stations, physical inputs and statistical thresholds
as R48fix1. Full input files are NOT committed; 202 input hashes passed preflight.
Formal D run reserves fusion weight 0.10 before searches. Partial shadow is off in
that run and is tested independently. The private FLOAT root is never conditioned
by a discarded search. Root snapshots and final-domain row/proof records are kept.

The selection benchmark uses the actual R48 00:30 history-conditioned domain,
not the unscreened original history pool. All gains and fixture risk parents are
held identical between new/reference algorithms. Only prefixes 64/128 are timed
with both algorithms; full 845 rows are timed with the new one. Compilation was
concurrent, so these numbers are diagnostic, not the formal GNSS speed ratio.

Scripts retain execution-machine absolute paths and require matching external
inputs/libraries. Freeze and formal startup remain gated on successful full build,
regression tests, actual PEA BLAS injection interception and a real epoch smoke.
