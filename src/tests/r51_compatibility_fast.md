# R51 exact integer compatibility in the existing product image

## Scope

Only the candidate loop in `r51SearchBlocks` changes. Its physical source store
is a local value, and its complete initial history has already been compiled
into `ZhangR51PhysicalImage`. All other physical compatibility call sites,
physical archive handling, tree/arc transport, LAMBDA search/limits, ratio=3,
numerical validity checks, posterior conditioning and writer authorization are
unchanged. Expensive diagnostic decompositions remain disabled.

## Exact equivalence

Let the initial integer domain be H*N=b and the searched physical targets be
T*N. The existing exact image proves

    {T*N : N integer, H*N=b} = {t0 + G*q : q integer}.

G spans the entire integer image (including nonprimitive divisibility), and the
rational left inverse L satisfies L*G=I. The current projector/offset represents
q=L*(T*N-t0). A reduced candidate B*q=z is lifted into physical rows by the
unchanged `zhangR51LiftRationalRows`, clearing denominators with nonzero integer
row scales D. On the initial domain, its residual is exactly D*(B*q-z).
Therefore a union of lifted physical candidates is feasible if and only if the
same union B*q=z is feasible over integer q. Hidden/retired physical variables
remain existential in the compiled image; none is set to zero or discarded.

The new per-call checker retains only accepted quotient rows. It assesses a
proposal read-only, then commits its ticket only after numerical conditioning
and decision-proof checks succeed. A rejection leaves its domain unchanged.
Tickets are bound to the checker and its generation, preventing stale/foreign
commits. The checker cannot be copied; it is destroyed at the end of one block
search and rebuilt for a different stage, epoch or physical image.

The old implementation rebuilt and decomposed the full physical union for
every candidate. The new implementation still performs exact integer
feasibility, but on the already available free-image coordinates. It avoids
both repeated dense physical matrices and their full physical-space right
transforms. It does not claim those checks can be removed globally.

## Verification

`zhang_r51_compatibility_tests` compares decisions against the unchanged
full-physical feasibility algorithm for deterministic and randomized domains:
hidden-arc parity, nonprimitive coupled images, rank-deficient targets, 130-bit
coefficients, inconsistent initial history, overlapping/contradictory new
integers, discarded-but-feasible trials and stale/foreign receipts. A bounded
synthetic benchmark reports the substep cost separately from domain compilation.

Production validation also runs the existing full-rank suite, actual production
LAMBDA/KF/transport contracts, and cross-day FLOAT/product/covariance replay plus
checkpoint restore. A prefix of the interrupted full-network run is retained
for comparison with the replacement run. Full-network timing must be reported
from actual completed timers, not extrapolated from a small benchmark.
