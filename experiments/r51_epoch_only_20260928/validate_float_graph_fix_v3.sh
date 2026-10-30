#!/usr/bin/env bash
set -euo pipefail

# Run only after building unit_tests: both production-object harnesses are
# dependencies of that target. Bind this fresh test log to the exact pea file
# subsequently frozen, and reject concurrent binary replacement.
repo=/mnt/c/Users/rx/Documents/GINAN/ginan-r51-ratio-only-nis-skip-20260923
build=/home/rx/GINAN/build-r51-handoff-20260930
[[ -x $repo/bin/pea && -x $build/graph_runtime_contract/graph_runtime_contract &&
   -x $build/float_qc_contract/float_qc_contract ]] || {
    echo 'build pea and unit_tests before validation' >&2; exit 2;
}
binary_before=$(sha256sum "$repo/bin/pea" | cut -d' ' -f1)
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
    ctest --test-dir "$build" --output-on-failure
binary_after=$(sha256sum "$repo/bin/pea" | cut -d' ' -f1)
[[ $binary_after == "$binary_before" ]] || {
    echo 'pea changed during validation; refusing success provenance' >&2; exit 2;
}
printf 'VALIDATED_PEA_SHA256=%s\n' "$binary_after"
