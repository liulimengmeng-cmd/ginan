FLOAT phase-scoped all-CPU profile (2026-10-19)

Default behaviour
  ZHANG_R51_FLOAT_THREADS unset or all: use all process-visible logical CPUs.
  Positive integer: controlled lower thread budget, at most available CPUs.
  legacy: retain inherited thread settings (does not undo matrix optimisations).
  Invalid values fail rather than silently selecting a different budget.

Scheduling
  Station preprocessing, per-receiver measurements and matrix-row assembly:
    OpenMP uses the selected budget; OpenBLAS and Eigen each use one thread.
  Joint network KF matrix products/solves:
    OpenBLAS uses the selected budget; Eigen's outer OpenMP pool is disabled.
  Innovation covariance factorisation:
    Cap BLAS workers to physical cores visible within the CPU affinity mask.
    This avoids the measured SMT panel-synchronisation penalty on Ryzen 9600X.
  Independent KF chunks retain the existing outer-parallel / BLAS-one rule.
  Guards restore pre-FLOAT effective OpenMP, Eigen and BLAS settings on all exits.
  No measurement model, covariance partition, ratio policy or NIS gate is changed.

Memory/computation changes
  Remove unused measurement-square identity.
  Build H_star D H_star' only for advanced_postfits, where it is actually needed.
  Materialise omega inverse/QinvH only in the existing omega-enabled branch.
  Compute P - KHP in its destination, avoiding a state-square KHP temporary.
  Symmetrise the posterior covariance once instead of two complete passes.
  Preserve vector measurement order; remove global Eigen reset-to-auto side effects.
  Reduce the measurement-assembly error flag safely across workers and never
  index a matrix with an undefined (-1) key.

Reproducible checks (Ubuntu-22.04-G)
  cmake --build /home/rx/GINAN/build-r51-handoff-20260930 --target pea -j2
  python3 src/tests/run_r51_float_qc.py /home/rx/GINAN/build-r51-handoff-20260930
  cmake --build /home/rx/GINAN/build-r51-handoff-20260930 \
    --target zhang_r51_float_parallel_contract_build -j2
  ctest --test-dir /home/rx/GINAN/build-r51-handoff-20260930 \
    -R zhang_r51_float_parallel_contract_tests --output-on-failure

The benchmark driver links the real PEA objects. Its independent small-matrix
oracle uses information-form estimation and checks the complete posterior for
eight combinations of Joseph, advanced-postfit and omega modes. Separate tests
check actual threads, no nested worker pools, exception restoration, measurement
order and correlated-noise assembly. The larger benchmark uses deterministic
synthetic inputs, not original station observations. It is NOT a real-network
epoch-speed result or scientific product-quality certification.

Paired benchmark audit
  python3 src/tests/audit_r51_float_parallel.py \
    --baseline /home/rx/GINAN/float_all_cpu_20261019/baseline/float_parallel_benchmark \
    --optimized /home/rx/GINAN/build-r51-handoff-20260930/float_parallel_contract/float_parallel_benchmark \
    --output /mnt/c/Users/rx/Documents/GINAN/artifacts/float_all_cpu_20261019/paired_final

The baseline executable was linked against pre-change production objects before
rebuilding PEA. Preserve it and all old products. The audit rotates run order,
runs sequentially, records CPU time/peak RSS and compares every posterior mean
and covariance entry. Do not run CPU-heavy jobs concurrently with timing.

No hour-long AR run is started by any of these tests. A future AR restart must
use the new compiled binary explicitly; old frozen launchers/binaries remain
unchanged. The original FLOAT checkpoint identity pins are not weakened. New
checkpoints record the FLOAT policy, setting, logical CPUs and physical cores.
