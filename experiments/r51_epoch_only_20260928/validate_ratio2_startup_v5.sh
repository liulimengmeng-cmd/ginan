#!/usr/bin/env bash
set -euo pipefail
repo=/mnt/c/Users/rx/Documents/GINAN/ginan-r51-ratio-only-nis-skip-20260923
evidence=/mnt/c/Users/rx/Documents/GINAN/artifacts/osb_ratio2_domain_budget_fix_20261001_v5
config_root=/home/rx/GINAN/frozen-r51-satellite-fix-20260918/config
[[ -x $repo/bin/pea && -d $evidence && ! -e $evidence/startup_ratio2_verbose.log &&
   ! -e $evidence/startup_invalid_ratio_verbose.log ]] || exit 2
configs=(
  "$config_root/zhang_global_2024199_180_base.yaml"
  "$config_root/zhang_global_2024199_180_inputs.yaml"
  "$config_root/zhang_global_2024199_180_product.yaml"
  "$config_root/zhang_global_2024199_180_e29_180network_product_1h.yaml"
  "$config_root/zhang_global_2024199_180_e29_service_30s_1h.yaml"
  "$config_root/case.yaml"
)
cd /mnt/d/GINAN_R20/inputData
env -u ZHANG_E29_RESTORE_BUNDLE ZHANG_R51_ENABLE=1 ZHANG_R51_RATIO_ONLY=1 \
  ZHANG_R51_RATIO_THRESHOLD=2 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  "$repo/bin/pea" --dry-run -y "${configs[@]}" -d case \
  > "$evidence/startup_ratio2_verbose.log" 2>&1
grep -q 'R51_RATIO_POLICY.*effective_threshold=2' "$evidence/startup_ratio2_verbose.log"
set +e
env -u ZHANG_E29_RESTORE_BUNDLE ZHANG_R51_ENABLE=1 ZHANG_R51_RATIO_ONLY=1 \
  ZHANG_R51_RATIO_THRESHOLD=invalid OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  "$repo/bin/pea" --dry-run -y "${configs[@]}" -d case \
  > "$evidence/startup_invalid_ratio_verbose.log" 2>&1
invalid_status=$?
set -e
[[ $invalid_status -ne 0 ]]
grep -q 'R51_RATIO_POLICY_INVALID' "$evidence/startup_invalid_ratio_verbose.log"
printf 'PASS production ratio2 startup and invalid-policy rejection; dry-run only\n'
