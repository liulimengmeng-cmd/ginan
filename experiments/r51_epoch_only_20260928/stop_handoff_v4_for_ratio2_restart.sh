#!/usr/bin/env bash
set -euo pipefail

# Stop this exact old experiment, never all PEA processes or the WSL guest.
target_pid=34
expected=/home/rx/GINAN/frozen-osb-handoff-budget-fix-20261001-v4/pea
old_run=/home/rx/GINAN/local_timing/r51_20261001_osb_branch_repair_ar1h_handoff_budget_fix_v4
[[ $(id -u) -eq 0 ]] || exit 2
[[ $(readlink "/proc/$target_pid/exe") == "$expected" ]] || {
    echo 'old PID identity changed; refusing to signal it' >&2; exit 2;
}
read -r -a initial_stat < "/proc/$target_pid/stat"
start_ticks=${initial_stat[21]}
printf 'OLD_PID=%s START_TICKS=%s\n' "$target_pid" "$start_ticks"
kill -TERM "$target_pid"
for attempt in {1..30}; do
    [[ -e /proc/$target_pid/exe ]] || break
    sleep 1
done
if [[ -e /proc/$target_pid/exe ]]; then
    [[ $(readlink "/proc/$target_pid/exe") == "$expected" ]] || exit 2
    read -r -a current_stat < "/proc/$target_pid/stat"
    [[ ${current_stat[21]} == "$start_ticks" ]] || exit 2
    echo 'Graceful stop timed out after 30 seconds; sending KILL to verified old PID'
    kill -KILL "$target_pid"
    for attempt in {1..10}; do
        [[ -e /proc/$target_pid/exe ]] || break
        sleep 1
    done
fi
[[ ! -e /proc/$target_pid/exe ]] || exit 2
printf 'OLD_PROCESS_STOPPED=1\n'
find "$old_run/outputs/zhang_internal_products.csv.epochs" -name COMMIT.json -printf '%h\n' | sort | tail -3
tail -20 "$old_run/runner.log"
