#!/usr/bin/env bash
set -euo pipefail
repo=/mnt/c/Users/rx/Documents/GINAN/ginan-r51-ratio-only-nis-skip-20260923
build=/home/rx/GINAN/build-r51-handoff-20260930
evidence=/mnt/c/Users/rx/Documents/GINAN/artifacts/osb_ratio2_domain_budget_fix_20261001_v5
[[ ! -e $evidence && ! -L $evidence ]] || { echo 'fresh evidence directory required' >&2; exit 2; }
mkdir -p "$evidence"
chmod +x "$repo/experiments/r51_epoch_only_20260928/cmake_windows_git/git"
export PATH="$repo/experiments/r51_epoch_only_20260928/cmake_windows_git:$PATH"
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 cmake -S "$repo/src" -B "$build" \
    -DGIT_EXECUTABLE="$repo/experiments/r51_epoch_only_20260928/cmake_windows_git/git" \
    2>&1 | tee "$evidence/configure.log"
# Only one compiler at a time: the live frozen AR process still owns most RAM.
nice -n 10 cmake --build "$build" --target pea -j 1 2>&1 | tee "$evidence/build_pea.log"
nice -n 10 cmake --build "$build" --target unit_tests -j 1 2>&1 | tee "$evidence/build_tests.log"
bash "$repo/experiments/r51_epoch_only_20260928/validate_float_graph_fix_v3.sh" \
    2>&1 | tee "$evidence/tests.log"
bash -n "$repo/experiments/r51_epoch_only_20260928/run_osb_user_strategy_restart.sh"
sha256sum "$repo/bin/pea" | tee "$evidence/binary.sha256"
printf 'BUILD_AND_CONTRACT_TESTS_COMPLETE; no AR process started or stopped\n'
