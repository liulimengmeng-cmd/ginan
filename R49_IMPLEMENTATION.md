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

A: implemented; production BLAS injection/replay and full tests pending build.

B1: full-set exact feasibility fast path with unchanged greedy fallback;
sparse incremental rational rank only after separate integer feasibility;
immutable proof-vector closure cache; zero-budget bridge returns before compile.
Original selector retained as a numerical/selection reference in tests.

B2: exact-posterior/order/epoch root-factor cache in full KF conditioner;
H-only decomposition reuse with fresh RHS/null checks; target marginal cache;
Gram/projection covariance retains full cross covariance. No tolerance changes.
