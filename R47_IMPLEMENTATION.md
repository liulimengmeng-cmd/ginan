# R47 implementation and acceptance record

Baseline: R46 source commit `05e5840`. Each stage is validated, committed and
pushed separately. Experimental coverage is not an acceptance substitute for
physical identity, integer feasibility, proof closure or final delivery checks.

## A: physical pullback

`zhangPullbackProductPhysicalRelations` projects the intended full physical
functional into the authoritative chart, then into the available posterior.
Both projections require exact re-expansion equality. The downstream chord
coordinates are rebound to the authoritative chart as well. Private subgraphs
need no square or unimodular inverse. Integers outside exact binary64 range
are rejected before numeric search. Existing affine constraint pullback is
retained, including `RT a = k - R d`.

Regression coverage: a pivot where the intended row is `[-1,1]` and the old
binding is `[0,1]`; a rectangular subgraph; missing posterior columns; invalid
physical support; exact numeric conversion. Existing full-state product gauge
and same-posterior integer/real gauge closure tests remain required.

## B: catalogue and available whole lattice

The runtime compiler retains all named expressions, affine structural identities,
integer recovery maps, and the saturated missing-column cancellation lattice.
Private trees use the design-approved shadow-only gate until their frontend
semantics are independently proven. The legacy named-star view remains separate;
stage E consumes the general search domain and publishes only proved consequences.

Actual design fixtures: 21:00 rowwise rank 16, whole rank 22; 27:30 catalogue
28, rank 22 with all six dependent identities preserved. Main executable builds;
411 test cases and 16,675 assertions pass. Structural zeros are not AR proofs.

## C: conditioning and immutable delivery

Product-image membership no longer filters state conditioners. The final
transaction freezes all admitted receipts, new equations, original decision
parents, exact physical expansions, ordered chart and phase-segment identities.
Delivery validates this immutable union on the protected FLOAT root and applies
the full-state equality conditioner, retaining cross-covariance updates. It no
longer performs a second Ledger selection. Infeasible integer affine systems
are rejected; nonprimitive but divisible equations are accepted.

415 tests / 16,700 assertions pass; main executable builds. Regression includes
the 1288/425 mixed-conditioner example, zero-cross-covariance control, full-state
updates, contradictory equalities and deletion of conditioner/proof dependencies.

## D: history selection and clean recertification source

Retired and unavailable physical integers are eliminated from the whole historical
lattice before strict physical projection. Surviving combinations retain every
contributing original proof. History subsets are selected before conditioning,
with deduplicated ancestor-union risk, exact affine feasibility and a reserved
quarter of the ceiling for new search. Redundant rows create no risk refund.
The recertification source is now loaded from the protected FLOAT KF owner, not
the potentially held-conditioned PRE_FRESH snapshot. Stage E executes both
controlled search paths using this actual root and their explicit parent sets.

Main executable and full regression pass, including retired-arc cancellation,
shared ancestors, reserved budget, and zero-increment-risk affine conflicts.

## E: official product quotient, route control and atomic delivery

The official joint product path now searches the exact affine integer image of
both complete available product lattices. It compacts only active columns,
preserves nonprimitive congruences, and uses global decorrelation followed by
sequential quotient blocks. Overlaps are retested as actual common functionals
from an uncontroversial pre-round posterior; joint integer feasibility and NIS
must pass before a round becomes an admitted conditioner.

FLOAT recertification and history-conditioned routes start from the same
protected FLOAT root. All attempted block allocations, including failed and
unchosen routes, charge the candidate-family budget. Immutable candidates freeze
the chosen route's conditioners, new constraints and original proof ancestry.
Named products require exact missing-column cancellation and exact consequence
proof; structural zeros alone never produce AR certificates.

The writer stages PRODUCT_FIXED covariance bytes and commits them only with the
physical Ledger transaction. Failed transactions discard fixed product rows and
restore both integer and gauge Ledger snapshots. Mixed integer evidence can
persist without pretending that it certifies a dual-frequency broadcast graph.
Diagnostics distinguish transported history, admitted conditioners, new integer
rank, product consequences, candidate graph, final published graph and ancestry
risk. Old pair records explicitly state retention/revocation checks.

Validation before final commit: main executable builds; 421 tests / 16,746
assertions pass. New regressions cover nonprimitive cosets, empty/infeasible
quotients, mixed conditioners, inactive network columns and sequential proof
receipts. A final post-commit source/build consistency check and full regression
are required before freezing the executable.

## Experiment and acceptance boundary

Start one fresh process only after A-E are pushed and the final build is verified.
Use 180 stations, 30-second sampling and the inclusive 2024-07-17 00:00-00:30
observation window (61 epochs), unchanged observation models and thresholds,
`product_relation_feedback: false`, and `transactional_integer_fixing: true`.
Disable the superseded sequential shadow route; E is the official route.

Actual 21:00 and 27:30 design fixtures are covered by exact regression. Full-run
14:00/19:00/21:00/26:00/27:30/29:00 behavior and comparison against R37/R45/R46
remain empirical acceptance items, not claims established by compilation or
unit tests. Private chart rebases remain shadow-only pending semantic proof.
Coverage, formal covariance and certified graph size do not independently prove
integer truth or positioning improvement. Preserve all raw outputs and logs.
