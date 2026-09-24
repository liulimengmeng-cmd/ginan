#!/usr/bin/env bash
set -euo pipefail

# Invoke only inside a private mount namespace as root. The original D: outputs
# remain visible and unchanged outside that namespace.
[[ $(id -u) -eq 0 ]] || { echo 'root mount namespace required' >&2; exit 2; }
[[ $(readlink /proc/self/ns/mnt) != $(readlink /proc/1/ns/mnt) ]] || {
  echo 'refusing to bind outputs in the host mount namespace' >&2
  exit 2
}
run_root=${ZHANG_LOCAL_TIMING_RUN_ROOT:-/home/rx/GINAN/local_timing/r51_20260924_instrumented_01}
run_timeout_seconds=${ZHANG_LOCAL_TIMING_TIMEOUT_SECONDS:-7200}
case "$run_timeout_seconds" in
  ''|*[!0-9]*) echo 'timeout must be a positive number of seconds' >&2; exit 2 ;;
esac
(( run_timeout_seconds > 0 )) || { echo 'timeout must be positive' >&2; exit 2; }
case "$run_root" in
  /home/rx/GINAN/local_timing/*) ;;
  *) echo 'run root must be inside /home/rx/GINAN/local_timing' >&2; exit 2 ;;
esac
configured_output=/mnt/d/GINAN_R20/inputData/outputs/zhang_r51_satfix_ratio3_ar1h_2024199_180_20260918
configured_checkpoints=/mnt/d/GINAN_R20/zhang_r51_satfix_ratio3_ar1h_2024199_180_20260918_checkpoints
binary=/home/rx/GINAN/build-r51-timing-20260924/bin/pea
config_root=/home/rx/GINAN/frozen-r51-satellite-fix-20260918/config
restore=/mnt/d/GINAN_R20/zhang_r51_ratio3_float2h_ar1h_2024199_180_20260917_checkpoints_source_before_ratioonly/E29-20240717015930-e240-rb06f1454a2e1

[[ -f "$binary" && -f "$restore/checkpoint.bundle" ]] || { echo 'binary or checkpoint missing' >&2; exit 2; }
[[ -d "$configured_output" && -d "$configured_checkpoints" ]] || { echo 'configured targets missing' >&2; exit 2; }
[[ $(findmnt -no FSTYPE -T /home/rx/GINAN) == ext4 ]] || { echo 'run root is not ext4' >&2; exit 2; }
install -d -o rx -g rx /home/rx/GINAN/local_timing
[[ ! -e "$run_root" ]] || { echo 'unique run directory already exists' >&2; exit 2; }
install -d -o rx -g rx "$run_root/outputs" "$run_root/checkpoints"
mount --bind "$run_root/outputs" "$configured_output"
mount --bind "$run_root/checkpoints" "$configured_checkpoints"
[[ $(stat -c %d:%i "$run_root/outputs") == $(stat -c %d:%i "$configured_output") ]]
[[ $(stat -c %d:%i "$run_root/checkpoints") == $(stat -c %d:%i "$configured_checkpoints") ]]

cd /mnt/d/GINAN_R20/inputData
exec runuser -u rx -- env \
  -u ZHANG_R51_EXPENSIVE_DIAGNOSTICS \
  -u ZHANG_R51_REUSE_RATIO_20240717 \
  -u ZHANG_R49_SNAPSHOT_DIRECTORY \
  -u ZHANG_R50_RATIO_ONLY \
  -u ZHANG_R49_PARTIAL_SHADOW \
  -u ZHANG_R49_BLAS_INJECT_INVALID \
  ZHANG_R51_ENABLE=1 \
  ZHANG_R51_RATIO_ONLY=1 \
  ZHANG_R51_REUSE_FLOAT_20240717=1 \
  ZHANG_R49_FUSION_WEIGHT=0.10 \
  ZHANG_P0_RESOURCE_PROBE=1 \
  OMP_NUM_THREADS=4 \
  OPENBLAS_NUM_THREADS=4 \
  MKL_NUM_THREADS=1 \
  ZHANG_R51_AR_START_GPST_SECONDS=1405216800 \
  ZHANG_E29_RESTORE_BUNDLE="$restore" \
  /usr/bin/time -v timeout --signal=TERM --kill-after=30s "${run_timeout_seconds}s" \
  "$binary" -q -y \
  "$config_root/zhang_global_2024199_180_base.yaml" \
  "$config_root/zhang_global_2024199_180_inputs.yaml" \
  "$config_root/zhang_global_2024199_180_product.yaml" \
  "$config_root/zhang_global_2024199_180_e29_180network_product_1h.yaml" \
  "$config_root/zhang_global_2024199_180_e29_service_30s_1h.yaml" \
  "$config_root/case.yaml" -d case > "$run_root/runner.log" 2>&1
