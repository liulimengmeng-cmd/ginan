# R50 ratio-only acceptance ablation

Baseline: R49 formal source 793814c422cceaada3e92d64f7175bc139c2c9d9.
Enable only with `ZHANG_R50_RATIO_ONLY=1`; default-off retains legacy acceptance.
Ratio = second-best weighted squared distance / best weighted squared distance, threshold 3.

- Integer estimation stays LAMBDA: same Z reduction, LDL moments, depth-first integer enumeration and nested suffix order. The bootstrap success cutoff no longer removes dimensions; failed ratio shrinks the same suffix. This may increase search time materially.
- Ratio requires two distinct returned integer vectors. The existing candidate-limit termination is preserved; no new claim of globally complete ILS enumeration is made. The ratio compares the best two returned candidates. Equal costs are retained in a multimap. Legacy enumeration is byte-preserved in the default-off branch. Missing second candidate rejects, equal-cost candidates reject, zero best with positive second accepts.
- NIS, bootstrap, marginal Perr and family-risk thresholds no longer veto active network/product/history decisions. Numerical validity, PSD/null-space consistency and exact affine integer compatibility still veto.
- New stochastic integer claims use the Ratio-tested LAMBDA paths. Historical confirmations retain integer identity, phase segments, epoch maturation and original proof dependencies; their additional NIS/Perr checks are diagnostic only.
- Search attempt limits, block ordering, bridge/fusion route allocation, exact tree pivots and integer transport remain unchanged. Nominal allocation constants still bound the number of attempts; accumulated historical risk no longer consumes their availability.
- All new R50 proof atoms explicitly carry the uninformative bound 1 and provenance `R50_RATIO_ONLY_NO_FFR_GUARANTEE`. A sum of these fields is diagnostic bookkeeping, not a probability claim. Raw NIS and bootstrap diagnostics remain visible. Product AR_VALID means this experimental policy only.
- Unchanged input window: 180 stations, 2024-07-17 00:00:00 through 00:30:00 GPST, 30-second intervals, 61 epochs. This is 30 minutes of observations, not a 30-minute wall-time timeout.
- Isolated worktree/config/output/checkpoint directories preserve R49. No input data copied to the experiment evidence bundle.

Validation includes top-two brute-force reference checks, default-off regression suite, R50 numerical/identity/risk policy tests, checkpoint tests, protected algorithm byte comparison and a one-epoch real PEA smoke run before formal launch.

Reference caveat: at a two-candidate test cap, one of 18 two-dimensional fixtures did not return the global top two under R49 candidate-limit termination. The formal R49/R50 configuration retains its candidate set size of 200; the reference harness also tests that setting. This known baseline search limitation is preserved intentionally to isolate the requested validation ablation; the experiment must not be presented as a globally complete ILS/Ratio benchmark.
