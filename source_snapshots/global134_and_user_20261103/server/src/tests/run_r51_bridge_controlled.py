"""Isolated authorized bridge validation run from the ORIGINAL pure-FLOAT checkpoint.

Run as WSL root: SCRIPT BINARY EVIDENCE_DIRECTORY SUFFIX [--stop-at HH:MM:SS]
Without --stop-at this is a full hour. A formal launch requires a successful
validation gate file; this script never converts a running test to a passed one.
"""
from pathlib import Path
from datetime import datetime,timedelta
import argparse
import hashlib
import json
import os
import pwd
import re
import shutil
import signal
import subprocess
import time

parser=argparse.ArgumentParser()
parser.add_argument('binary', type=Path)
parser.add_argument('evidence', type=Path)
parser.add_argument('suffix')
parser.add_argument('--stop-at')
parser.add_argument('--formal-gates',type=Path)
args=parser.parse_args()
assert os.geteuid()==0, 'Private mount namespace needs WSL root'
assert re.fullmatch('[A-Za-z0-9_-]+',args.suffix)
assert bool(args.stop_at) != bool(args.formal_gates), 'Select controlled stop or accepted formal gates'
if args.stop_at:
    assert re.fullmatch('02:[0-5][0-9]:[03]0',args.stop_at)
else:
    gates=json.loads(args.formal_gates.read_text())
    assert gates.get('all_required_tests_passed') is True
    assert gates.get('event_0206_verified') is True and gates.get('event_0237_verified') is True

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()

source=args.binary.resolve(strict=True)
assert source.name=='pea' and os.access(source,os.X_OK)
evidence=args.evidence.resolve()
assert evidence.parent==Path('/home/rx/GINAN/bridge_continuity_20261020')
assert not evidence.exists()
assert subprocess.run(['pgrep','-x','pea'],stdout=subprocess.DEVNULL).returncode!=0,'Another solver exists'
checkpoint=Path('/mnt/d/GINAN_R20/zhang_r51_ratio3_float2h_ar1h_2024199_180_20260917_checkpoints_source_before_ratioonly/E29-20240717015930-e240-rb06f1454a2e1/checkpoint.bundle')
assert digest(checkpoint)=='ab397e655774da5c0a9632becd08fe171508b93ffeaf2dd7b04de899bae2182a'
checkpoint_manifest=checkpoint.parent/'checkpoint_manifest.json'
checkpoint_manifest_hash=digest(checkpoint_manifest)
config=Path('/home/rx/GINAN/frozen-cycle-slip-20261018-qualified-gf-pinned-v1/config')
config_target=Path('/home/rx/GINAN/frozen-r51-satellite-fix-20260918/config')
launcher=Path('/home/rx/GINAN/frozen-osb-ratio2-domain-budget-fix-20261001-v5/launcher.sh')
run=Path('/home/rx/GINAN/local_timing')/('r51_20261001_osb_branch_repair_ar1h_'+args.suffix)
assert not run.exists() and launcher.is_file() and config.is_dir()
evidence.mkdir()
binary=evidence/'pea'
shutil.copy2(source,binary)
binary_hash=digest(binary)
if args.formal_gates:
    assert gates['binary_sha256']==binary_hash,'Gates belong to a different binary'
account=pwd.getpwnam('rx')
capture=evidence/'bridge_inputs'
capture.mkdir()
os.chown(capture,account.pw_uid,account.pw_gid)
events=[]
state=dict(stage='controlled_event_validation' if args.stop_at else 'formal_ar_hour',
    status='STARTING',source_binary=str(source),binary=str(binary),binary_sha256=binary_hash,run_root=str(run),
    checkpoint=str(checkpoint),checkpoint_sha256=digest(checkpoint),ratio=2,joint_nis=False,
    checkpoint_manifest_sha256=checkpoint_manifest_hash,launcher_sha256=digest(launcher),
    epoch_only=True,stop_at=args.stop_at,config_sha256={p.name:digest(p) for p in config.glob('*.yaml')},
    events=events,required_events=['02:06:00','02:37:30'],scientific_acceptance=False)
def save():
    temporary=evidence/'run_state.tmp'
    temporary.write_text(json.dumps(state,indent=2))
    temporary.replace(evidence/'run_state.json')
save()
started=time.monotonic()
environment={**{k:v for k,v in os.environ.items() if not k.startswith('ZHANG_')},
    'ZHANG_R51_EVIDENCE_SLIP_QC':'1','ZHANG_R51_SLIP_RAW_WARMUP':'1',
    'ZHANG_R51_BRIDGE_REPLAY_DIRECTORY':str(capture),
    'ZHANG_R51_BRIDGE_CAPTURE_EPOCHS':'02:03:30,02:06:00,02:27:00,02:37:00,02:37:30,02:54:30,02:56:30',
    'ZHANG_R51_BLOCK_CONDITIONER':'incremental','ZHANG_R51_FLOAT_MATRIX_MODE':'symmetric',
    'ZHANG_R51_HANDOFF_FLOAT_RETENTION_SECONDS':'300','ZHANG_R51_HANDOFF_EVENT_BUDGET_MS':'60000'}
state['pinned_environment']={key:value for key,value in environment.items() if key.startswith('ZHANG_') and
    key not in ('ZHANG_R51_BRIDGE_REPLAY_DIRECTORY',)}
state['handoff_float_retention_seconds']=environment.get('ZHANG_R51_HANDOFF_FLOAT_RETENTION_SECONDS','300')
state['bridge_event_budget_ms']=environment.get('ZHANG_R51_HANDOFF_EVENT_BUDGET_MS','60000')
if args.formal_gates:
    for key in ('config_sha256','pinned_environment','checkpoint_sha256','checkpoint_manifest_sha256','launcher_sha256'):
        assert gates[key]==state[key],'Formal run identity differs from validated candidate: '+key
with (evidence/'launch.log').open('w') as output:
    process=subprocess.Popen(['unshare','--mount','--propagation','private','--','bash','-c',
        'set -euo pipefail; ulimit -c 0; mount --bind "$1" "$2"; '
        'mount --bind "$6" "$6"; mount -o remount,bind,ro "$6"; exec bash "$3" "$4" "$5"',
        'bridge-run',str(config),str(config_target),str(launcher),str(binary),str(run),str(checkpoint.parent)],
        stdout=output,stderr=subprocess.STDOUT,env=environment)
    state['launcher_pid']=process.pid;save()
    solver=None;interrupted=False
    def children(pid):
        result=[pid]
        for item in result:
            try:result.extend(map(int,Path(f'/proc/{item}/task/{item}/children').read_text().split()))
            except FileNotFoundError:pass
        return result
    def owned(pid):
        try:return Path(f'/proc/{pid}/exe').resolve()==binary
        except OSError:return False
    try:
        while process.poll() is None:
            if solver is None:
                candidates=[pid for pid in children(process.pid) if owned(pid)]
                assert len(candidates)<=1
                if candidates:
                    solver=candidates[0]
                    assert digest(Path(f'/proc/{solver}/exe'))==binary_hash
                    state.update(solver_pid=solver,status='RUNNING');save()
                    print(json.dumps(dict(event='verified_solver',pid=solver,run_root=str(run))),flush=True)
            committed=[]
            for path in (run/'outputs/zhang_internal_products.csv.epochs').glob('*/COMMIT.json'):
                try:
                    record=json.loads(path.read_text())
                    if record.get('ar') is True:committed.append((record['epoch'],str(path)))
                except (OSError,json.JSONDecodeError):pass
            committed.sort()
            if committed:
                state.update(committed_ar_epochs=len(committed),last_committed_epoch=committed[-1][0])
                if args.stop_at and committed[-1][0][-8:]>=args.stop_at and not interrupted:
                    assert solver is not None and owned(solver)
                    os.kill(solver,signal.SIGINT);interrupted=True
                    events.append(dict(event='controlled_stop_after_commit',epoch=committed[-1][0]))
            state['elapsed_seconds']=time.monotonic()-started;save()
            if state['elapsed_seconds']>36*3600:
                raise TimeoutError('Authorized run exceeded 36-hour safety limit')
            time.sleep(5)
        result=process.wait()
    except BaseException:
        if solver is not None and owned(solver):os.kill(solver,signal.SIGKILL)
        if process.poll() is None:process.terminate()
        state.update(status='FAILED_OR_INTERRUPTED',scientific_acceptance=False);save()
        raise
state.update(returncode=result,status='CONTROLLED_STOPPED' if interrupted and result==130 else
    'EXITED' if result==0 else 'FAILED',solver_still_running=solver is not None and owned(solver))
final_committed=[]
for path in (run/'outputs/zhang_internal_products.csv.epochs').glob('*/COMMIT.json'):
    record=json.loads(path.read_text())
    if record.get('ar') is True:final_committed.append(record['epoch'])
final_committed.sort()
state.update(committed_ar_epochs=len(final_committed),
    last_committed_epoch=final_committed[-1] if final_committed else None)
save()
print(json.dumps(state,indent=2),flush=True)
assert digest(checkpoint)==state['checkpoint_sha256'] and digest(checkpoint_manifest)==checkpoint_manifest_hash, \
    'Original checkpoint changed during isolated run'
assert result==(130 if args.stop_at else 0)
assert not state['solver_still_running']
if not args.stop_at:
    expected=[(datetime(2024,7,17,2)+timedelta(seconds=30*index)).strftime('%Y-%m-%d %H:%M:%S')
        for index in range(120)]
    assert final_committed==expected, \
        'Formal run exited without all 120 committed AR epochs'
