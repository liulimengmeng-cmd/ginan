#!/usr/bin/env bash
set -euo pipefail
repo=/mnt/c/Users/rx/Documents/GINAN/ginan-r51-ratio-only-nis-skip-20260923
build=/home/rx/GINAN/build-r51-handoff-20260930
evidence=/mnt/c/Users/rx/Documents/GINAN/artifacts/osb_handoff_budget_fix_20261001
mkdir -p "$evidence"
cmake -S "$repo/src" -B "$build" -DGIT_EXECUTABLE=/mnt/d/Git/cmd/git.exe \
    2>&1 | tee "$evidence/configure_v4.log"
if [[ ${1:-} == budget-only ]]; then
    cmake --build "$build" --target zhang_r51_handoff_search_budget_tests \
        zhang_r51_current_bridge_representation_tests -j 2
    "$repo/bin/zhang_r51_handoff_search_budget_tests"
    "$repo/bin/zhang_r51_current_bridge_representation_tests"
    exit
fi
cmake --build "$build" --target pea -j 2 2>&1 | tee "$evidence/build_pea_v4.log"
cmake --build "$build" --target unit_tests -j 2 2>&1 | tee "$evidence/build_tests_v4.log"
bash "$repo/experiments/r51_epoch_only_20260928/validate_float_graph_fix_v3.sh" \
    2>&1 | tee "$evidence/tests_final_v4.log"
