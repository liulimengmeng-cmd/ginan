#!/usr/bin/env bash
set -euo pipefail

# Launch as root in a private mount namespace. The first argument is the newly
# built pea binary; the second is a fresh run directory on the WSL ext4 volume.
# Keep the full one-hour case configuration. The caller may stop monitoring
# after two committed epochs while this process continues in the background.
if [[ $# -ne 2 ]]; then
  echo "usage: $0 /absolute/path/to/pea /home/rx/GINAN/local_timing/r51_20261001_osb_branch_repair_ar1h_ID" >&2
  exit 2
fi

binary=$1
run_root=$2
configured_output=/mnt/d/GINAN_R20/inputData/outputs/zhang_r51_satfix_ratio3_ar1h_2024199_180_20260918
configured_checkpoints=/mnt/d/GINAN_R20/zhang_r51_satfix_ratio3_ar1h_2024199_180_20260918_checkpoints
restore=/mnt/d/GINAN_R20/zhang_r51_ratio3_float2h_ar1h_2024199_180_20260917_checkpoints_source_before_ratioonly/E29-20240717015930-e240-rb06f1454a2e1
config_root=/home/rx/GINAN/frozen-r51-satellite-fix-20260918/config

[[ $(id -u) -eq 0 ]] || { echo 'root required' >&2; exit 2; }
[[ $(readlink /proc/self/ns/mnt) != $(readlink /proc/1/ns/mnt) ]] || {
  echo 'private mount namespace required' >&2; exit 2;
}
[[ $binary == /* && -x $binary ]] || { echo 'an absolute executable pea path is required' >&2; exit 2; }
[[ $run_root == /home/rx/GINAN/local_timing/r51_20261001_osb_branch_repair_ar1h_* ]] || {
  echo 'run directory must use the dedicated local_timing prefix' >&2; exit 2;
}
[[ $(readlink -f -- "$(dirname -- "$run_root")") == /home/rx/GINAN/local_timing ]] || {
  echo 'run directory parent must be /home/rx/GINAN/local_timing' >&2; exit 2;
}
[[ ! -e $run_root && ! -L $run_root ]] || { echo 'run directory already exists' >&2; exit 2; }
[[ -f $restore/checkpoint.bundle && -f $restore/checkpoint_manifest.json ]] || {
  echo 'restore checkpoint is incomplete' >&2; exit 2;
}
[[ -d $configured_output && -d $configured_checkpoints ]] || {
  echo 'configured mount targets are missing' >&2; exit 2;
}
[[ $(findmnt -no FSTYPE -T /home/rx/GINAN) == ext4 ]] || {
  echo 'run root must be on ext4' >&2; exit 2;
}
grep -Eq '^[[:space:]]*end_epoch: 2024-07-17 02:59:30$' "$config_root/case.yaml" || {
  echo 'the full one-hour end epoch is missing from case.yaml' >&2; exit 2;
}
grep -Eq '^[[:space:]]*product_mode: HOU_OSB_LIKE$' "$config_root/case.yaml" &&
grep -Eq '^[[:space:]]*hou_product_coordinate: HYBRID_STABLE$' "$config_root/case.yaml" || {
  echo 'the final case config must select HOU_OSB_LIKE/HYBRID_STABLE for single-basis handoff' >&2
  exit 2
}

install -d -o rx -g rx "$run_root" "$run_root/outputs" "$run_root/checkpoints"
mount --bind "$run_root/outputs" "$configured_output"
mount --bind "$run_root/checkpoints" "$configured_checkpoints"
[[ $(stat -c %d:%i "$run_root/outputs") == $(stat -c %d:%i "$configured_output") ]] || {
  echo 'output bind mount verification failed' >&2; exit 2;
}
[[ $(stat -c %d:%i "$run_root/checkpoints") == $(stat -c %d:%i "$configured_checkpoints") ]] || {
  echo 'checkpoint bind mount verification failed' >&2; exit 2;
}

printf 'RUN_ROOT=%s\nBINARY=%s\n' "$run_root" "$binary"
sha256sum "$binary"
cd /mnt/d/GINAN_R20/inputData
exec runuser -u rx -- env \
  -u ZHANG_R51_EXPENSIVE_DIAGNOSTICS \
  -u ZHANG_R51_REUSE_RATIO_20240717 \
  -u ZHANG_R51_REUSE_FLOAT_THREADS_8_1 \
  -u ZHANG_R51_AR_THREADS \
  -u ZHANG_R49_SNAPSHOT_DIRECTORY \
  -u ZHANG_R50_RATIO_ONLY \
  -u ZHANG_R49_PARTIAL_SHADOW \
  -u ZHANG_R49_BLAS_INJECT_INVALID \
  ZHANG_R51_ENABLE=1 \
  ZHANG_R51_RATIO_ONLY=1 \
  ZHANG_R51_RATIO_THRESHOLD=2 \
  ZHANG_R51_EPOCH_ONLY_AR=1 \
  ZHANG_R51_SINGLE_BASIS_HANDOFF=1 \
  ZHANG_R51_PRODUCT_REPAIR=1 \
  ZHANG_R51_HANDOFF_SEARCH_BUDGET_MS=10000 \
  ZHANG_R51_INTERN_HANDOFF_PROOFS=1 \
  ZHANG_R51_REUSE_FLOAT_20240717=1 \
  ZHANG_R49_FUSION_WEIGHT=0.10 \
  ZHANG_P0_RESOURCE_PROBE=1 \
  OMP_NUM_THREADS=4 \
  OPENBLAS_NUM_THREADS=4 \
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
