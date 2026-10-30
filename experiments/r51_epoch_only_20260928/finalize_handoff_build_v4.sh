#!/usr/bin/env bash
set -euo pipefail
repo=/mnt/c/Users/rx/Documents/GINAN/ginan-r51-ratio-only-nis-skip-20260923
build=/home/rx/GINAN/build-r51-handoff-20260930
evidence=/mnt/c/Users/rx/Documents/GINAN/artifacts/osb_handoff_budget_fix_20261001
chmod +x "$repo/experiments/r51_epoch_only_20260928/cmake_windows_git/git"
export PATH="$repo/experiments/r51_epoch_only_20260928/cmake_windows_git:$PATH"
cmake -S "$repo/src" -B "$build" \
    -DGIT_EXECUTABLE="$repo/experiments/r51_epoch_only_20260928/cmake_windows_git/git" \
    2>&1 | tee "$evidence/configure_v4_versioned.log"
cmake --build "$build" --target pea -j 2 2>&1 | tee "$evidence/build_pea_v4_versioned.log"
cmake --build "$build" --target unit_tests -j 2 2>&1 | tee "$evidence/build_tests_v4_versioned.log"
if [[ -f $evidence/tests_final_v4.log ]]; then
    [[ ! -e $evidence/tests_before_metadata_repair_v4.log ]] || exit 2
    cp "$evidence/tests_final_v4.log" "$evidence/tests_before_metadata_repair_v4.log"
fi
bash "$repo/experiments/r51_epoch_only_20260928/validate_float_graph_fix_v3.sh" \
    2>&1 | tee "$evidence/tests_final_v4.log"
