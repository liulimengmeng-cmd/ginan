# R51 ratio-only certificate entailment speed fix

The third AR epoch spent hours in `ZhangR51PhysicalStore::entails`, rebuilding a physical affine domain for every gauge certificate. The gauge scan does not mutate its physical archive or current chart until all candidate families have been collected.

The scan now owns one immutable copy of H and b. It lazily computes the same exact integer particular solution and saturated kernel once, then checks each target with t*K=0 and t*x0=rhs. The old direct-row receipt shortcut is preserved. Unknown nonzero physical variables fail entailment. No process/global cache is retained: new epochs, arcs, right-hand sides and call contexts create new snapshots. No integer search, ratio threshold, NIS policy, transport rule or precision tolerance is changed.

`ZHANG_R51_GAUGE_ENTAILMENT_CACHE` reports queries, direct queries, affine builds and build time. New WL/L1 phase timers separate ILS, physical compatibility, joint HNF, root NIS and covariance conditioning; these add observation only.

Validation: `zhang_r51_entailment_tests` compares 739 queries with the previous production algorithm, including nonprimitive constraints, inconsistent lattices, hidden variables, 130-bit coefficients, deterministic random systems, changed RHS and arc retirement. A 64-variable/96-query fixture verified 96 repeated reference queries against one cached affine build; observed 0.449 s vs 0.0041 s on this host. This benchmark is not a whole-network speed forecast.

Other AR costs remain: per-block full compatibility/HNF/NIS and full KF factorization are unchanged. The first epoch's uninstrumented 35-minute region requires the new timers before choosing the next optimization.
