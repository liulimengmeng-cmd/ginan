# R49 implementation and acceptance

Source: OSB / 分析解算速度并测试, conversation 6aa68697-6fac-83e9-be18-acaf20169431,
latest design 9e70982f-4190-4cda-9a00-42d9194d40d0. Full four-turn transcript retained locally.
Baseline 820c2cb; preserve R48 execution snapshot 94af7df and all old outputs.

Stages: A observable search exits/source provenance/actual PEA BLAS audit;
B1 identical-order incremental exact selection and budget-first bridge;
B2 posterior/conditioner/target cache and square-root marginals;
C dynamic cuts, exact target equivalence, atomic append, whole-route finalization;
D prebudgeted target recertification from alternate route;
E partial bridge shadow only, never public dual AR authorization.

Acceptance: unchanged statistical thresholds; integer infeasibility and missing-column
cancellation regressions; rejected append preserves prior domain; unselected route
targets get explicit outcomes; local NIS passing does not override whole NIS.
Formal run: same 180 stations, 2024-07-17 00:00:00--00:30:00, 61 epochs,
OMP=4 / BLAS=1, unique output/checkpoint roots, source and binary frozen.
Do not begin formal run until compilation, numerical/contract tests and startup pass.

A: implemented; the final production binary intercepts the invalid-DGEMV injection
with exit 86 and a call stack. Runtime OpenBLAS reports LP64.

B1: full-set exact feasibility fast path with unchanged greedy fallback;
sparse incremental rational rank only after separate integer feasibility;
immutable proof-vector closure cache; zero-budget bridge returns before compile.
Original selector retained as a numerical/selection reference in tests.

B2: exact-posterior/order/epoch root-factor cache in full KF conditioner;
H-only decomposition reuse with fresh RHS/null checks; target marginal cache;
Gram/projection covariance retains full cross covariance. No tolerance changes.

C: condition-domain metadata, exact target-on-domain compiler and signed/translated
canonical requests; dynamic component reconstruction after commit; attempted-risk
retained, trial append atomically replaces only after full route gates. Initial
routes are finalized before comparison. Joint gate reports base/increment/whole.
The sparse rank basis never authorizes integer membership. Physical target maps
still use the existing missing-column cancellation and final candidate writer.

D: fusion search share is reserved before either route (default 0 for control;
formal fusion ablation ZHANG_R49_FUSION_WEIGHT=0.10). Maximum 4 recertifications.
Only target identities transfer, never alternate RHS, moments or certificates.
All alternate accepted bridge targets and pair targets get explicit outcomes.
Every accepted addition rebuilds the final full domain, risk and product catalogue.

E: optional ZHANG_R49_PARTIAL_SHADOW counterfactual on rejected rank-2 bridges;
maps PAR rows/RHS back through exact projector and offsets, audits full joint NIS,
classifies WL/L1/mixed/full consequences. Never mutates formal domain or Ledger.
Shadow timing and risk are separate; disabled for formal performance comparison.

Audit completion: unique posterior serials; one optional binary network-root
snapshot (mean, column-major covariance, exact column order) per root; final
domain records include every ordinary-search and bridge row/RHS, source and
full statistical/algebraic dependency lists. Formal run enables root snapshots.
Actual PEA invalid-DGEMV injection returned 86 with argument values and call stack.
OpenBLAS runtime ABI is checked against USE64BITINT when that backend is linked.
Reference selector tests: 320 randomized cases + integer parity regressions pass.
Recorded R48 00:30 subset fixture: 128-row selections identical; full 845 rows
accepted by new selector. Timing is diagnostic under concurrent compilation.


Validation: 433/433 full regression cases and 9/9 checkpoint cases passed.
All three standalone reference suites passed with EIGEN_USE_BLAS=1, the actual
OpenBLAS archive and the DGEMV audit wrapper. A real 180-station first-epoch
smoke is a separate mandatory gate before freezing and launching the formal run.
These checks do not establish a full-window speedup or improved GNSS quality.

Real startup caught DGEMV M=0, N=2588, LDA=0 in final product-domain NIS.
The final NIS call now checks an empty domain before Eigen evaluates H*x;
empty domains receive no statistical acceptance or certificate. The regression
reproduces the 0x2588 shape under the BLAS audit and compares nonempty NIS.
The failed smoke log is retained locally; a fresh output root is required for retry.
