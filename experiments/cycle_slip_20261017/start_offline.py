"""Dispatch a frozen, test-qualified offline run; never stop another solver."""
import argparse
from datetime import datetime, timedelta, timezone
import json
import os
from pathlib import Path
import subprocess
import time


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--frozen',type=Path,required=True)
    args=parser.parse_args()
    if os.geteuid()!=0:raise SystemExit('Launch dispatcher must run as root')
    frozen=args.frozen.resolve()
    if frozen.parent!=Path('/home/rx/GINAN') or not frozen.name.startswith('frozen-cycle-slip-'):
        raise SystemExit('Unexpected frozen directory')
    manifest=json.loads((frozen/'manifest.json').read_text())
    if manifest['restart_authorization']!='USER_REQUEST_AFTER_TESTS_PASSED':raise SystemExit('No restart authorization')
    if not json.loads((frozen/'validation_result.json').read_text())['all_execution_checks_passed']:
        raise SystemExit('Readiness gate failed')
    if (frozen/'launch_status.json').exists():raise SystemExit('Launch already dispatched; inspect it, do not duplicate')
    if Path(manifest['run_root']).exists():raise SystemExit('Run output already exists')
    other=subprocess.run(['pgrep','-x','pea'],capture_output=True,text=True)
    if other.returncode==0:raise SystemExit('Another PEA exists; do not duplicate or terminate it: '+other.stdout)
    with (frozen/'launch.log').open('w') as log:
        process=subprocess.Popen(['unshare','--mount','--propagation','private','bash',
            str(frozen/'launch_offline.sh'),str(frozen)],stdin=subprocess.DEVNULL,
            stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    time.sleep(1)
    status={'status':'LAUNCH_DISPATCHED' if process.poll() is None else 'LAUNCH_FAILED',
        'dispatcher_pid':process.pid,'exit_code':process.poll(),
        'wall_time_china':datetime.now(timezone(timedelta(hours=8))).isoformat(),
        'run_root':manifest['run_root'],'binary_sha256':manifest['binary_sha256'],
        'note':'Dispatch is not proof of a restored checkpoint, committed epoch, or product quality'}
    (frozen/'launch_status.json').write_text(json.dumps(status,indent=2))
    print(json.dumps(status))
    return 0 if process.poll() is None else 1


if __name__=='__main__':raise SystemExit(main())
