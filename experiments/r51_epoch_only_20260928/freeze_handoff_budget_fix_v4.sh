#!/usr/bin/env bash
set -euo pipefail

# Freeze only. The caller keeps a foreground WSL exec alive when launching;
# printing a detached PID alone was not reliable on this host in v3.
repo=/mnt/c/Users/rx/Documents/GINAN/ginan-r51-ratio-only-nis-skip-20260923
frozen=/home/rx/GINAN/frozen-osb-handoff-budget-fix-20261001-v4
run_root=/home/rx/GINAN/local_timing/r51_20261001_osb_branch_repair_ar1h_handoff_budget_fix_v4
launcher=experiments/r51_epoch_only_20260928/run_osb_user_strategy_restart.sh
validated_log=/mnt/c/Users/rx/Documents/GINAN/artifacts/osb_handoff_budget_fix_20261001/tests_final_v4.log

[[ $# -eq 1 && $1 =~ ^[a-f0-9]{64}$ ]] || {
    echo 'supply the SHA256 of the production binary whose tests passed' >&2; exit 2;
}
validated_binary_sha=$1
[[ $(id -u) -eq 0 ]] || { echo 'run as root' >&2; exit 2; }
[[ ! -e $frozen && ! -L $frozen && ! -e $run_root && ! -L $run_root ]] || {
    echo 'v4 destination exists; refusing overwrite' >&2; exit 2;
}
[[ -x $repo/bin/pea && -f $repo/$launcher ]] || exit 2
[[ $(sha256sum "$repo/bin/pea" | cut -d' ' -f1) == "$validated_binary_sha" ]] || {
    echo 'production binary changed after validation' >&2; exit 2;
}
grep -Eq '^100% tests passed, 0 tests failed out of [0-9]+$' "$validated_log"
tr -d '\r' < "$validated_log" | grep -Fxq "VALIDATED_PEA_SHA256=$validated_binary_sha"

install -d -o rx -g rx "$frozen"
install -m 755 -o rx -g rx "$repo/bin/pea" "$frozen/pea"
[[ $(sha256sum "$frozen/pea" | cut -d' ' -f1) == "$validated_binary_sha" ]] || exit 2
install -m 755 -o rx -g rx "$repo/$launcher" "$frozen/launcher.sh"
tar -czf "$frozen/source.tar.gz" -C "$repo" src "$launcher" \
    experiments/r51_epoch_only_20260928/freeze_handoff_budget_fix_v4.sh \
    experiments/r51_epoch_only_20260928/build_handoff_budget_fix_v4.sh \
    experiments/r51_epoch_only_20260928/finalize_handoff_build_v4.sh \
    experiments/r51_epoch_only_20260928/cmake_windows_git/git \
    experiments/r51_epoch_only_20260928/validate_float_graph_fix_v3.sh
chown rx:rx "$frozen/source.tar.gz"
bash -n "$frozen/launcher.sh"
sha256sum "$frozen/pea" "$frozen/launcher.sh" "$frozen/source.tar.gz"
printf 'FROZEN=%s\nRUN_ROOT=%s\nREADY_TO_LAUNCH=1\n' "$frozen" "$run_root"
