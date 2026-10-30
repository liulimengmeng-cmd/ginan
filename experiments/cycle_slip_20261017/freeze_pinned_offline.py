"""Freeze tested QC with the original, exactly pinned target configuration.

The optional receiver TRACE addition in freeze_offline changes canonical target
identity. Do not relax the C++ migration whitelist for extra logging: restore
the complete original config directory and explicitly record reduced tracing.
This wrapper changes no detector, tested binary, source fixture or checkpoint.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


def sha(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--validation', type=Path, required=True)
    parser.add_argument('--frozen', type=Path, required=True)
    parser.add_argument('--run-root', type=Path, required=True)
    parser.add_argument('--restart-after-passed-tests', action='store_true')
    args = parser.parse_args()
    if not args.restart_after_passed_tests:
        raise SystemExit('Restart after passed tests must be requested')
    original = Path('/home/rx/GINAN/frozen-r51-satellite-fix-20260918/config')
    original_hashes = {p.name: sha(p) for p in sorted(original.glob('*.yaml'))}
    if not original_hashes or 'case.yaml' not in original_hashes:
        raise SystemExit('Pinned target config is incomplete')
    subprocess.run([sys.executable, str(Path(__file__).with_name('freeze_offline.py')),
        '--validation', str(args.validation), '--frozen', str(args.frozen),
        '--run-root', str(args.run_root), '--restart-after-passed-tests'], check=True)
    manifest_path = args.frozen / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    if manifest['status'] != 'FROZEN_NOT_STARTED':
        raise SystemExit('Unexpected frozen status; do not change a running experiment')
    if original_hashes != {p.name: sha(p) for p in sorted(original.glob('*.yaml'))}:
        raise SystemExit('Original config changed while freezing')
    # These are new, not-yet-launched artifact copies, never original inputs.
    for name in original_hashes:
        shutil.copy2(original / name, args.frozen / 'config' / name)
    frozen_hashes = {p.name: sha(p) for p in sorted((args.frozen / 'config').glob('*.yaml'))}
    if frozen_hashes != original_hashes:
        raise SystemExit('Pinned config copy mismatch; do not launch')
    manifest['config_hashes'] = frozen_hashes
    manifest['trace_change'] = 'NONE: exact original pinned config; no added receiver TRACE'
    manifest['freeze_wrapper_sha256'] = sha(Path(__file__))
    manifest['freeze_wrapper'] = str(Path(__file__).resolve())
    manifest['config_policy'] = 'EXACT_ORIGINAL_TARGET_BYTES; CPP_MIGRATION_WHITELIST_UNCHANGED'
    manifest['limitations'].append(
        'Per-link QC TRACE is not newly enabled in this network run; detailed evidence remains in the validated raw-reader replay logs')
    manifest_path.write_text(json.dumps(manifest, indent=2))
    print(json.dumps({'frozen': str(args.frozen), 'run_root': str(args.run_root),
        'config_policy': manifest['config_policy'], 'status': manifest['status']}))


if __name__ == '__main__':
    main()
