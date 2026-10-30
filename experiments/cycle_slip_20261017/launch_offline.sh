#!/usr/bin/env bash
set -euo pipefail
# Run only inside a new private mount namespace, as root. The old launcher and
# configuration remain unchanged outside this namespace.
[[ $# -eq 1 ]] || { echo 'usage: launch_offline.sh /absolute/frozen-directory' >&2; exit 2; }
frozen=$1
[[ $frozen == /home/rx/GINAN/frozen-cycle-slip-* ]] || { echo 'unexpected frozen directory' >&2; exit 2; }
[[ $(id -u) -eq 0 && $(readlink /proc/self/ns/mnt) != $(readlink /proc/1/ns/mnt) ]] || {
    echo 'root in a private mount namespace required' >&2; exit 2;
}
run_root=$(python3 - "$frozen" <<'PY'
import hashlib,json,sys
from pathlib import Path
root=Path(sys.argv[1]);data=json.loads((root/'manifest.json').read_text())
def sha(path):
    value=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):value.update(block)
    return value.hexdigest()
assert data['status']=='FROZEN_NOT_STARTED' and data['restart_authorization']=='USER_REQUEST_AFTER_TESTS_PASSED'
assert json.loads((root/'validation_result.json').read_text())['all_execution_checks_passed']
assert sha(root/'pea')==data['binary_sha256']
assert sha(Path(data['checkpoint']))==data['checkpoint_sha256']
assert sha(Path(data['base_launcher']))==data['base_launcher_sha256']
for path,digest in data['input_hashes'].items():assert sha(Path(path))==digest
for name,digest in data['config_hashes'].items():assert sha(root/'config'/name)==digest
print(data['run_root'])
PY
)
export ZHANG_R51_EVIDENCE_SLIP_QC=1
export ZHANG_R51_SLIP_RAW_WARMUP=1
mount --bind "$frozen/config" /home/rx/GINAN/frozen-r51-satellite-fix-20260918/config
exec bash /home/rx/GINAN/frozen-osb-ratio2-domain-budget-fix-20261001-v5/launcher.sh "$frozen/pea" "$run_root"
