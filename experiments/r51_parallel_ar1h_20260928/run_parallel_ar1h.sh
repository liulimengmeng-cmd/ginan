#!/usr/bin/env bash
set -euo pipefail

# Launch with: wsl -d Ubuntu-22.04-G -u root -- unshare --mount --propagation private bash <this-file>
# Bind only the frozen configuration's output paths inside this mount namespace.
[[ $(id -u) -eq 0 ]] || { echo 'root required' >&2; exit 2; }
[[ $(readlink /proc/self/ns/mnt) != $(readlink /proc/1/ns/mnt) ]] || {
  echo 'private mount namespace required' >&2; exit 2;
}

run_root=/home/rx/GINAN/local_timing/r51_20260928_parallel_ar1h_03
configured_output=/mnt/d/GINAN_R20/inputData/outputs/zhang_r51_satfix_ratio3_ar1h_2024199_180_20260918
configured_checkpoints=/mnt/d/GINAN_R20/zhang_r51_satfix_ratio3_ar1h_2024199_180_20260918_checkpoints
restore=/mnt/d/GINAN_R20/zhang_r51_ratio3_float2h_ar1h_2024199_180_20260917_checkpoints_source_before_ratioonly/E29-20240717015930-e240-rb06f1454a2e1
binary=/home/rx/GINAN/build-r51-p0p4-20260927/bin/pea
config_root=/home/rx/GINAN/frozen-r51-satellite-fix-20260918/config

[[ ! -e "$run_root" ]] || { echo 'run directory already exists' >&2; exit 2; }
[[ -f "$binary" && -f "$restore/checkpoint.bundle" ]] || {
  echo 'binary or checkpoint missing' >&2; exit 2;
}
[[ -d "$configured_output" && -d "$configured_checkpoints" ]] || {
  echo 'configured targets missing' >&2; exit 2;
}
[[ $(findmnt -no FSTYPE -T /home/rx/GINAN) == ext4 ]] || {
  echo 'run root must be ext4' >&2; exit 2;
}

install -d -o rx -g rx "$run_root/outputs" "$run_root/checkpoints"
mount --bind "$run_root/outputs" "$configured_output"
mount --bind "$run_root/checkpoints" "$configured_checkpoints"
[[ $(stat -c %d:%i "$run_root/outputs") == $(stat -c %d:%i "$configured_output") ]]
[[ $(stat -c %d:%i "$run_root/checkpoints") == $(stat -c %d:%i "$configured_checkpoints") ]]

cd /mnt/d/GINAN_R20/inputData
exec runuser -u rx -- env \
  -u ZHANG_R51_EXPENSIVE_DIAGNOSTICS \
  -u ZHANG_R51_REUSE_RATIO_20240717 \
  -u ZHANG_R51_EPOCH_ONLY_AR \
  -u ZHANG_R49_SNAPSHOT_DIRECTORY \
  -u ZHANG_R50_RATIO_ONLY \
  -u ZHANG_R49_PARTIAL_SHADOW \
  -u ZHANG_R49_BLAS_INJECT_INVALID \
  ZHANG_R51_ENABLE=1 \
  ZHANG_R51_RATIO_ONLY=1 \
  ZHANG_R51_REUSE_FLOAT_20240717=1 \
  ZHANG_R51_REUSE_FLOAT_THREADS_8_1=1 \
  ZHANG_R49_FUSION_WEIGHT=0.10 \
  ZHANG_P0_RESOURCE_PROBE=1 \
  ZHANG_R51_AR_THREADS=8 \
  OMP_NUM_THREADS=8 \
  OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 \
  ZHANG_R51_AR_START_GPST_SECONDS=1405216800 \
  ZHANG_E29_RESTORE_BUNDLE="$restore" \
  /usr/bin/time -v timeout --signal=TERM --kill-after=30s 14d \
  "$binary" -q -y \
  "$config_root/zhang_global_2024199_180_base.yaml" \
  "$config_root/zhang_global_2024199_180_inputs.yaml" \
  "$config_root/zhang_global_2024199_180_product.yaml" \
  "$config_root/zhang_global_2024199_180_e29_180network_product_1h.yaml" \
  "$config_root/zhang_global_2024199_180_e29_service_30s_1h.yaml" \
  "$config_root/case.yaml" -d case > "$run_root/runner.log" 2>&1
