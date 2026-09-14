#!/usr/bin/env bash
set -Eeuo pipefail
audit=/mnt/c/Users/rx/Documents/GINAN/r50_work_20260914
data=/mnt/d/GINAN_R20/inputData
frozen=/home/rx/GINAN/frozen-r50-ratio-only-20260914
case_name=zhang_p3_exp3_G3_r50_ratio_only_2024199_180_003000_20260914
output="$data/outputs/$case_name"
log="$audit/$case_name.runner.log"
state="$audit/$case_name.state.txt"
exec 8>/home/rx/GINAN/.ginan_global_pea_experiment.lock
flock -n 8 || exit 88
pgrep -x pea >/dev/null && exit 90
[[ ! -e /mnt/d/GINAN_R20/r50_root_snapshots_20260914 && ! -e "$output" && ! -e "$log" && ! -e "$state" && ! -e /mnt/d/GINAN_R20/r50_checkpoints_20260914 ]] || exit 89
cd "$frozen"
sha256sum -c SHA256SUMS >"$audit/r50_preflight_hashes.log" || exit 93
root_source=$(findmnt -nro SOURCE /)
root_options=$(findmnt -nro OPTIONS /)
[[ ",$root_options," == *,rw,* ]] || exit 92
mem_available_kib=$(awk '/MemAvailable:/ {print $2}' /proc/meminfo)
[[ "$mem_available_kib" -ge 11500000 ]] || exit 94
record() { printf '%s %s\n' "$(date --iso-8601=seconds)" "$*" | tee -a "$state" "$log"; }
trap 'rc=$?; record "RUNNER_EXIT status=$rc"' EXIT
python3 "$audit/verify_r50.py" --preflight >"$audit/r50_root_preflight.json"
python3 - "$audit" <<'PY'
import json,sys
from pathlib import Path
w=Path(sys.argv[1])
p=json.loads((w/'r50_input_preflight.json').read_text())
assert p['files']==202 and not p['failures']
t=json.loads((w/'frozen_r50.json').read_text())
assert t['tests_passed'] and t['source_matches_build']
PY
record "PREFLIGHT PASS root=$root_source mem_available_kib=$mem_available_kib frozen_hashes=PASS inputs=202/202 tests=PASS expected_epochs=61"
record "R50_POLICY ratio_only=1 ratio_threshold=3 bootstrap_gate=0 nis_gate=0 perr_gate=0 risk_budget_gate=0 risk_bound_valid=0 exact_integer_checks=1 tree_transport_unchanged=1"
record "START case=$case_name binary_sha256=$(sha256sum "$frozen/bin/pea" | awk '{print $1}') OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 observation_window=2024-07-17T00:00:00/00:30:00 official_product_quotient=1 authoritative_feedback=0"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export ZHANG_R50_RATIO_ONLY=1
export ZHANG_R49_FUSION_WEIGHT=0.10
export ZHANG_R49_SNAPSHOT_DIRECTORY=/mnt/d/GINAN_R20/r50_root_snapshots_20260914
unset ZHANG_R49_PARTIAL_SHADOW ZHANG_R49_BLAS_INJECT_INVALID
cd "$data"
set +e
/usr/bin/time -v "$frozen/bin/pea" -q -y \
    "$frozen/config/zhang_global_2024199_180_base.yaml" \
    "$frozen/config/zhang_global_2024199_180_inputs.yaml" \
    "$frozen/config/zhang_global_2024199_180_product.yaml" \
    "$frozen/config/zhang_global_2024199_180_e29_180network_product_1h.yaml" \
    "$frozen/config/zhang_global_2024199_180_e29_service_30s_1h.yaml" \
    "$frozen/config/case.yaml" -d "$case_name" >>"$log" 2>&1
pea_status=$?
set -e
record "FINISH pea_status=$pea_status output_preserved=1"
set +e
python3 "$audit/verify_r50.py" --pea-status "$pea_status" >"$audit/r50_final_verification.json"
verification_status=$?
set -e
record "FINAL_VERIFICATION status=$verification_status report=$audit/r50_final_verification.json"
exit "$verification_status"
