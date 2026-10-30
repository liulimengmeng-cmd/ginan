#!/usr/bin/env bash
set -euo pipefail

repo=/mnt/c/Users/rx/Documents/GINAN/ginan-r51-ratio-only-nis-skip-20260923
frozen=/home/rx/GINAN/frozen-osb-ratio2-domain-budget-fix-20261001-v5
run_root=/home/rx/GINAN/local_timing/r51_20261001_osb_branch_repair_ar1h_ratio2_domain_budget_fix_v5
launcher=experiments/r51_epoch_only_20260928/run_osb_user_strategy_restart.sh
validation=/mnt/c/Users/rx/Documents/GINAN/artifacts/osb_ratio2_domain_budget_fix_20261001_v5
expected_sha=72f32aea698358a39a093613947b46a6808f22a302d774ef811f68ff3ec51be3

[[ $(id -u) -eq 0 ]] || { echo 'root required' >&2; exit 2; }
[[ ! -e $frozen && ! -L $frozen && ! -e $run_root && ! -L $run_root ]] || {
    echo 'v5 destination already exists; refusing overwrite' >&2; exit 2;
}
[[ $(sha256sum "$repo/bin/pea" | cut -d' ' -f1) == "$expected_sha" ]] || {
    echo 'validated binary changed' >&2; exit 2;
}
tr -d '\r' < "$validation/tests.log" | grep -Fxq '100% tests passed, 0 tests failed out of 33'
tr -d '\r' < "$validation/tests.log" | grep -Fxq "VALIDATED_PEA_SHA256=$expected_sha"
bash -n "$repo/$launcher"

install -d -o rx -g rx "$frozen"
install -m 755 -o rx -g rx "$repo/bin/pea" "$frozen/pea"
install -m 755 -o rx -g rx "$repo/$launcher" "$frozen/launcher.sh"
install -m 644 -o rx -g rx "$validation/validation_record.json" "$frozen/validation_record.json"
install -m 644 -o rx -g rx "$validation/tests.log" "$frozen/tests.log"
[[ $(sha256sum "$frozen/pea" | cut -d' ' -f1) == "$expected_sha" ]] || exit 2
tar -czf "$frozen/source.tar.gz" -C "$repo" src "$launcher" \
    experiments/r51_epoch_only_20260928/freeze_ratio2_domain_budget_v5.sh \
    experiments/r51_epoch_only_20260928/build_ratio2_domain_budget_v5.sh \
    experiments/r51_epoch_only_20260928/finish_ratio2_domain_budget_v5.sh \
    experiments/r51_epoch_only_20260928/validate_ratio2_startup_v5.sh \
    experiments/r51_epoch_only_20260928/cmake_windows_git/git \
    experiments/r51_epoch_only_20260928/validate_float_graph_fix_v3.sh
chown rx:rx "$frozen/source.tar.gz"
sha256sum "$frozen/pea" "$frozen/launcher.sh" "$frozen/source.tar.gz"
printf 'FROZEN=%s\nRUN_ROOT=%s\nREADY_TO_LAUNCH=1\n' "$frozen" "$run_root"
