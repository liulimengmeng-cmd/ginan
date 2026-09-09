# Experiment 5: CODE phase-OSB nonzero prior sensitivity

## Frozen design

Baseline: `codex/code-products-pppar`, commit
`91c1ceae9ecf1192419d1db83841d1751ef7cb37`, worktree
`D:/tec/ginan-main-code-products`. No estimator source changes are required.
The frozen PEA SHA-256 is
`3ab177d7ee7c06936d02384368fe0f090ccee0a16d4e6b6049fadde8b486b23c`.

Run root: `D:/tec/code-exp5-sigma-20260909-v2`.
The first preparation directory, `D:/tec/code-exp5-sigma-20260909`, is retained:
preparation reached fixture generation but failed when Linux Git encountered the
Windows-native linked-worktree path. No PEA run was launched there. The driver now
translates that metadata path; it never changes the worktree's `.git` file.

The predeclared grid is 0.003 and 0.010 m, each with a full 2880-epoch AR run and
a paired FLOAT probe-only run. Existing experiment 4c/4d is the zero-sigma control.
Neither success threshold 0.9999 nor ratio threshold 3 is changed. The six stations,
GPS L1C/L2W, G01 exclusion, 15-degree mask, receiver satellite-single-difference
integer coordinates, product means, 2024-07-17 observations, and 14:00-24:00 analysis
window are unchanged. RTS remains enabled, but coordinate comparison uses forward
`AR` and `PPP` blocks, as in the baseline.

## Model and narrow implementation

For each GPS satellite and signal, a single bias state has prior mean equal to the
original CODE OSB and prior variance sigma squared. These priors are independent
across satellite/signal keys. The state is shared across all six receivers and
persists through the day with zero process noise and no receiver chunking.
No cross-signal or cross-satellite covariance is inferred from OSB mean values.
The sigma choices are sensitivity hypotheses, not calibrated CODE product accuracy.

`scripts/experiment5_code_osb_sigma.py` makes explicitly named experimental BIA
fixtures. It replaces only bytes 92:103 (zero-based, standard-deviation field)
of the 64 GPS satellite L1C/L2W OSB lines. All other bytes, including means,
intervals, code OSBs, and headers, remain identical. The accompanying fixture
receipts explicitly state that these are not published CODE products. Original
CODE files are never edited. G01's fixture entries exist, but G01 remains excluded,
so the expected number of used bias states is 62.

Sigma is encoded in ns to seven decimal places. The receipt records the resulting
effective metre sigma and hash rather than assuming exact decimal representation.
The existing parser converts ns to metres, squares it, and the existing persistent
state path uses that covariance exactly once. No variance is repeatedly injected
as phase observation noise. Code-bias uncertainty treatment stays unchanged.

Each configuration is a standalone copy with exactly one BIA path. This avoids
GINAN's append behavior for input lists. `configuration_control_audit.json`
verifies that only the include path and selected BIA path differ from 4c/4d.
The original configuration comments/metadata are retained for bytewise comparison;
the enclosing experiment-5 run labels and plan identify the actual sigma condition.

## Source anchors and interpretation corrections

- `src/cpp/common/biasSINEXread.cpp:162`: ns-to-metre conversion;
  lines 187-192: standard deviation read and squared.
- `src/cpp/pea/ppp_obs.cpp:1420`: persistent external-bias condition;
  lines 1436-1442: product variance overrides configured initial P;
  lines 1452-1456: once-only state creation marker.
- `src/cpp/pea/ppp_ambres.cpp:1158`: analytic innovation covariance and NIS;
  line 1248: actual `filterKalman` call; lines 1293-1297: full-state discrepancy norm.
- `src/cpp/common/algebra.cpp:3060`: actual posterior screening;
  `src/cpp/pea/ppp_callbacks.cpp:66`: measurement-noise row/column deweighting.

The prior report's `1.399 m` shadow discrepancy label is incorrect. That diagnostic
is the norm of the whole state vector, combining different physical units, and
cannot be interpreted as a position displacement. Actual filtering also performs
posterior checks and may alter the update through callbacks; the analytic shadow
does not apply those callbacks. The supplemental audit therefore separates calls
with and without screening failures. It does not call this difference a proven
nonlinear-model error or attribute it solely to zero OSB uncertainty.

## Required audit outputs

For every run: return code, command/environment receipt, 2880 epochs, candidate and
failure counts, submitted rank, station/network feedback incidence, NIS distribution,
shadow norms conditional on screening, exact-once bias keys, all 62 state creations,
and effective prior covariance. FLOAT must contain zero feedback submissions.

For each pair: exact epoch coverage and 14:00-24:00 UTC 3D/horizontal/vertical RMS
against the same weekly IGS coordinate SINEX, including individual stations. No
station velocity propagation is inferred. No tuning to coordinates or alternate
burn-in window is performed. A lower NIS median alone does not certify consistency:
tail behavior, candidate selection, repeated constraints, and missing integer truth
remain relevant. Feedback, candidate, `SUBMITTED_UNVERIFIED`, verified fixing, and
coordinate improvement are separate quantities.

Every original input is enumerated with size and SHA-256 in `all_inputs.json`;
the final report includes the complete data table. The 23 files are verified both
before and after execution. Raw traces and the original experimental directories
are retained in place.

## Reproduction

Run under WSL `Ubuntu-22.04-G` from the code-products worktree, using a new output
root (the driver refuses to replace existing roots or run directories):

```sh
python3 scripts/experiment5_code_osb_sigma.py prepare /mnt/d/tec/NEW_EXP5_ROOT
python3 scripts/experiment5_code_osb_sigma.py run /mnt/d/tec/NEW_EXP5_ROOT --label sigma_3mm_ar --epochs 5
python3 scripts/experiment5_code_osb_sigma.py run /mnt/d/tec/NEW_EXP5_ROOT --label sigma_3mm_ar
python3 scripts/experiment5_code_osb_sigma.py run /mnt/d/tec/NEW_EXP5_ROOT --label sigma_3mm_float
python3 scripts/experiment5_code_osb_sigma.py run /mnt/d/tec/NEW_EXP5_ROOT --label sigma_10mm_ar
python3 scripts/experiment5_code_osb_sigma.py run /mnt/d/tec/NEW_EXP5_ROOT --label sigma_10mm_float
python3 scripts/audit_experiment5_code_osb_sigma.py /mnt/d/tec/NEW_EXP5_ROOT
```

The driver sets Boost's library path and OMP/OpenBLAS/MKL thread counts to one.
Four independent processes were run concurrently with approximately 12 GiB
available before launch; each used about 0.7 GiB resident memory during filtering.
No pre-existing PEA process was present or interrupted.
