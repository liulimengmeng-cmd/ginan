#!/usr/bin/env bash
set -euo pipefail
repo=/mnt/c/Users/rx/Documents/GINAN/ginan-r51-ratio-only-nis-skip-20260923
build=/home/rx/GINAN/build-r51-handoff-20260930
evidence=/mnt/c/Users/rx/Documents/GINAN/artifacts/osb_ratio2_domain_budget_fix_20261001_v5
[[ -d $evidence && ! -e $evidence/build_tests_retry.log && ! -e $evidence/tests.log ]] || exit 2
export PATH="$repo/experiments/r51_epoch_only_20260928/cmake_windows_git:$PATH"
nice -n 10 cmake --build "$build" --target unit_tests -j 1 \
    2>&1 | tee "$evidence/build_tests_retry.log"
bash "$repo/experiments/r51_epoch_only_20260928/validate_float_graph_fix_v3.sh" \
    2>&1 | tee "$evidence/tests.log"
bash -n "$repo/experiments/r51_epoch_only_20260928/run_osb_user_strategy_restart.sh"
sha256sum "$repo/bin/pea" | tee "$evidence/binary.sha256"
printf 'BUILD_AND_CONTRACT_TESTS_COMPLETE; no AR process started or stopped\n'
