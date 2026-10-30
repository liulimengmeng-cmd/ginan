#!/usr/bin/env bash
set -euo pipefail

# Execute only after the new production build and regression suite pass.
# Preserve all previous frozen binaries, runs, checkpoints and raw inputs.
repo=/mnt/c/Users/rx/Documents/GINAN/ginan-r51-ratio-only-nis-skip-20260923
frozen=/home/rx/GINAN/frozen-osb-float-graph-fix-20261001-v3
run_root=/home/rx/GINAN/local_timing/r51_20261001_osb_branch_repair_ar1h_float_graph_fix_v3
launcher=experiments/r51_epoch_only_20260928/run_osb_user_strategy_restart.sh
validated_log=/mnt/c/Users/rx/Documents/GINAN/artifacts/osb_float_graph_fix_20261001/tests_final_v3.log

[[ $# -eq 1 && $1 =~ ^[a-f0-9]{64}$ ]] || {
    echo 'supply the SHA256 of the production binary whose tests just passed' >&2; exit 2;
}
validated_binary_sha=$1

[[ $(id -u) -eq 0 ]] || { echo 'run as root for the private mount namespace' >&2; exit 2; }
[[ ! -e $frozen && ! -L $frozen && ! -e $run_root && ! -L $run_root ]] || {
    echo 'v3 frozen/run destination already exists; refusing overwrite' >&2; exit 2;
}
[[ -x $repo/bin/pea && -f $repo/$launcher ]] || {
    echo 'built pea or audited launcher is missing' >&2; exit 2;
}
[[ $(sha256sum "$repo/bin/pea" | cut -d' ' -f1) == "$validated_binary_sha" ]] || {
    echo 'built binary changed after validation; refusing launch' >&2; exit 2;
}
grep -Fq '100% tests passed, 0 tests failed out of 27' "$validated_log" || {
    echo 'the current v3 complete-suite success log is missing' >&2; exit 2;
}
tr -d '\r' < "$validated_log" | grep -Fxq "VALIDATED_PEA_SHA256=$validated_binary_sha" || {
    echo 'test success log is not bound to this production binary' >&2; exit 2;
}

install -d -o rx -g rx "$frozen"
install -m 755 -o rx -g rx "$repo/bin/pea" "$frozen/pea"
[[ $(sha256sum "$frozen/pea" | cut -d' ' -f1) == "$validated_binary_sha" ]] || {
    echo 'frozen binary does not match tested binary; refusing launch' >&2; exit 2;
}
install -m 755 -o rx -g rx "$repo/$launcher" "$frozen/launcher.sh"
tar -czf "$frozen/source.tar.gz" -C "$repo" src "$launcher" \
    experiments/r51_epoch_only_20260928/freeze_and_start_float_graph_fix_v3.sh \
    experiments/r51_epoch_only_20260928/validate_float_graph_fix_v3.sh
chown rx:rx "$frozen/source.tar.gz"
bash -n "$frozen/launcher.sh"
sha256sum "$frozen/pea" "$frozen/launcher.sh" "$frozen/source.tar.gz"

# Persist beyond the launching shell; no monitoring job or new schedule is
# created. The launcher maintains separate ext4 output/checkpoint mounts.
nohup unshare -m --propagation private bash "$frozen/launcher.sh" \
    "$frozen/pea" "$run_root" > "$frozen/launcher.log" 2>&1 < /dev/null &
launch_pid=$!
printf 'LAUNCHER_PID=%s\nFROZEN=%s\nRUN_ROOT=%s\n' "$launch_pid" "$frozen" "$run_root"
