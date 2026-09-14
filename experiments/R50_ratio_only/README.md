# R50 reproducibility records

Formal entry point: `/mnt/c/Users/rx/Documents/GINAN/r50_work_20260914/run_r50_30m.sh`.
The runner requires the immutable frozen binary/source bundle and refuses existing result directories. It runs the same 61 epochs and performs one final audit on exit; it does not poll or schedule monitoring.

The five base input/model YAML files are identical to R49. `case.yaml` changes only run identifiers and checkpoint locations. `ZHANG_R50_RATIO_ONLY=1` in the runner is the experimental treatment. Input data are referenced by hash and never copied here.

See repository `R50_IMPLEMENTATION.md` for acceptance semantics and the preserved baseline candidate-search limitation. Never interpret AR_VALID from this run as a 0.001 failure probability guarantee.
