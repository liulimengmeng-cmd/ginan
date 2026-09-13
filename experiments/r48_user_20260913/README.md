# R48 independent user PPP-AR experiment

The actual experiment uses the frozen 2024-07-17 00:03:00--00:24:30 PRODUCT_FIXED prefix (44 epochs). MARS, DYNG, NICO, BREW and JPLM are absent from the exact 180-station network input. Each user estimates position freely with identical products and full covariance; only the AR mode and feedback flag differ within each pair. Acceptance thresholds remain unchanged.

## Execution and provenance

Work: `/mnt/c/Users/rx/Documents/GINAN/r48_user_work_20260913_prefix`.
Binary: the unchanged frozen R48fix1 binary, source 94af7df, SHA256 9b70f939569a72564f332f8ea22d1c9adc83b38d9be38c19a948ca1cd5e415d0.

Preparation sequence: `prepare.py`, `enable_if_capture.py`, `restore_if_design.py`. Preparation refuses overwrites; an intentional rerun requires a new work/output root. Final case specification: `experiment_if.json`.

Run `run.py --spec experiment_if.json --case CASE` for individual cases, or omit `--case` for serial execution. Existing successful receipts are validated and skipped; a failed case stops the batch and is never silently overwritten. BREW's failure was retained, and JPLM was subsequently run via individual case selection. All user processes use one OMP/BLAS thread and niceness 10. The server is untouched.

Run `analyse.py` after all ten case attempts. It requires contiguous actual epochs, verifies POS references against IGS SINEX, checks formal product hashes, retains failures and exports separate WL/L1/committed-constraint evidence. `diagnose_brew.py` reproduces the failed BREW transaction with informational logs enabled in a unique output root.

## Results and limits

Four station pairs completed; BREW FLOAT and AR both failed at 00:06:30. Every case had zero actually committed user integer rank. Three stations admitted partial WL targets but none passed conditional L1. No positioning benefit was demonstrated. The full server R48 experiment continues independently; this user trial consumes only the completed prefix.

Raw input data are referenced and hashed, never copied. ENU is relative to the IGS weekly SINEX apriori (also used for initialization with a 100 m position prior), not an external ambiguity truth. The upstream five DGEMV warnings are unresolved. Failed setup attempts and the detailed BREW reproduction are retained in the work and output directories.

The `results` and `executed_config` folders preserve the final evidence and exact used configurations without input data or frozen product copies.
