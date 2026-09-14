# R51 implementation and experiment protocol

Source: OSB / 分析解算速度并测试, assistant message 9426892a-9739-4059-94d2-712c7f7a6f80. The latest two design turns were read and saved in R51_GPT_DESIGN.md. Linked attachments were not retrieved; no claim is made that their reference test archive was executed.

Base: R49 formal commit 793814c422cceaada3e92d64f7175bc139c2c9d9. Isolated branch codex/r51-certified-blocks. Activation requires ZHANG_R51_ENABLE=1; normal R49 behavior remains the default.

## Implemented behavior

- A: one physical compatibility view retains old integer variables; proposals carry exact rational-conflict or divisibility witnesses. Product ledger preflight operates on a copy and returns a privately constructed immutable receipt binding the numerical root, physical epoch, full ledger rows, RHS, metadata and proof DAG. New product constraints undergo a physical union check before candidate freeze. Current-epoch history-only product proposals are checked for exact entailment and zero new decision atoms. Failed enhancement restores product runtime and ledgers, attempts that current-root baseline, and otherwise emits current FLOAT products. New network HELD rows are persisted only after successful product commit.
- B: healthy historical blocks pass one whole-block NIS monitor and reuse original integers/proofs directly; exceptional failing blocks retain the existing isolation path. The 7.5e-4 historical reserve restriction becomes the existing 1e-3 certificate scope. Pending row confirmation uses stored RHS after a whole-block monitor. Historical risk is retained, never replaced with zero.
- C: block search uses the exact physical affine integer image, including existential old variables and nonprimitive congruences. A rational left inverse is audited exactly; constraints are lifted by exact denominator clearing. WL is followed by L1 even when no new WL succeeds. Initial blocks 32, at most 8 proposals per block, 128 actual ILS calls per stage, 2 scans; positive progress does not stop subsequent blocks. Covariance is the marginal after accepted conditions. Complete top-two ILS enumeration has a fixed 2,000,000-node limit; incomplete results are rejected.
- D: the final product search base is frozen after network L1 and includes all accepted network conditions. Existing R49 product bridge/route code is retained. Nonidentity tree transport is permitted only after an exact physical cycle round trip in both directions. Historical gauge mappings additionally require physical entailment; a matching satellite/segment label alone grants no integer authority.
- Publication: per-epoch private directories contain products, covariance, original proof metadata, COMMIT.json and SHA256SUMS. Renaming the directory commits the bundle. Aggregate CSV files are derived convenience views and do not provide cross-file atomicity. Physical events precede product snapshots and are not rolled back by product enhancement failure.
- Checkpoint schema 3 retains the full R51 physical archive and immutable proof DAG; old schema payloads must not silently restore as R51 history. This does not establish cross-version replay equivalence.

## Formal experiment

R49_MATCHED is the selected policy. Original success/ratio and whole-domain NIS checks remain enabled. Block thresholds may tighten to fit remaining certificate risk; they are not relaxed to produce fixing. R51_CONDITIONAL is an optional later ablation, not mixed into this run. The bounded adaptive block procedure is not presented as an empirically calibrated family-wise failure guarantee.

Window: 2024-07-17 00:00:00 through 01:00:00 GPST, inclusive, 121 epochs, 180 stations. OMP=4, OpenBLAS=1, MKL=1. Input data remains in its existing location. All output, checkpoint, root snapshot and frozen executable paths are new R51 paths.

Before launch: final-source compilation, legacy and R51 regression runs, exact-witness tests, physical-image tests, independent ILS oracle, incomplete-search rejection, input hash audit and actual PEA smoke. The running executable is a frozen copy and must never be overwritten.

## Evaluation boundaries

The formal run must specifically examine 00:25:30 transport and 00:29:30/00:30:00 conflicts. Unit tests and startup success do not prove those historical failures are resolved. Keep strict dual-frequency AR satellites, connected components, fixed network rank, actual ILS calls, history fallback, physical conflicts and exclusive phase timings separate. AR_VALID is internal certification, not independent integer truth or demonstrated PPP-AR positioning improvement.
