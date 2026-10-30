"""Authorized serial workflow: controlled validation, audit, then formal AR.

The formal hour is never launched after an exception, nonzero test exit, missing
gate, user interruption, or changed candidate. This is a one-off job, not a
monitoring schedule. Run as WSL root after build/CTest and candidate freeze.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import signal
import subprocess
import sys
import time

parser=argparse.ArgumentParser()
parser.add_argument('candidate',type=Path)
parser.add_argument('build',type=Path)
parser.add_argument('identity')
args=parser.parse_args()
candidate=args.candidate.resolve(strict=True)
base=candidate.parent
assert base==Path('/home/rx/GINAN/bridge_continuity_20261020')
assert os.geteuid()==0
identity=args.identity
assert identity and all(c.isalnum() or c in '-_' for c in identity)
control=base/('controlled_'+identity)
formal=base/('formal_'+identity)
state_path=base/('workflow_'+identity+'.json')
assert not state_path.exists() and not control.exists() and not formal.exists()
scripts=candidate/'source/src/tests'
binary=candidate/'pea'
manifest=json.loads((candidate/'manifest.json').read_text())
assert hashlib.sha256(binary.read_bytes()).hexdigest()==manifest['binary_sha256']
state=dict(status='STARTING',pid=os.getpid(),candidate=str(candidate),binary_sha256=manifest['binary_sha256'],
    controlled_evidence=str(control),formal_evidence=str(formal),formal_ar_started=False,
    ratio=2,joint_nis=False,epoch_only=True,started_unix=time.time())

def save():
    temporary=state_path.with_suffix('.tmp')
    temporary.write_text(json.dumps(state,indent=2))
    temporary.replace(state_path)

def execute(stage,arguments):
    state.update(status=stage,stage_started_unix=time.time());save()
    with (base/('workflow_'+identity+'_'+stage+'.log')).open('w') as stream:
        process=subprocess.Popen([sys.executable,*map(str,arguments)],stdout=stream,stderr=subprocess.STDOUT)
        state['child_pid']=process.pid;save()
        try:
            while process.poll() is None:
                if stage=='FORMAL_HOUR' and (formal/'run_state.json').exists():
                    live=json.loads((formal/'run_state.json').read_text())
                    state['formal_ar_started']='solver_pid' in live
                    state['formal_solver_pid']=live.get('solver_pid')
                    state['formal_last_committed_epoch']=live.get('last_committed_epoch')
                    save()
                time.sleep(5)
            code=process.wait()
        except BaseException:
            if process.poll() is None:
                process.send_signal(signal.SIGINT)
                process.wait(timeout=30)
            raise
    state['last_returncode']=code;save()
    if code:raise RuntimeError(stage+' failed with exit '+str(code))

def stop_requested(signum,frame):raise KeyboardInterrupt('Workflow stop signal '+str(signum))
signal.signal(signal.SIGTERM,stop_requested)
save()
try:
    execute('CONTROLLED_EVENT_REPLAY',[scripts/'run_r51_bridge_controlled.py',binary,control,
        'bridge_'+identity+'_controlled','--stop-at','02:38:00'])
    execute('EVENT_AND_ACCEPTANCE_AUDIT',[scripts/'audit_r51_bridge_controlled.py',control,args.build])
    gates=json.loads((control/'formal_gates.json').read_text())
    assert gates['all_required_tests_passed'] is True and gates['binary_sha256']==state['binary_sha256']
    execute('FORMAL_HOUR',[scripts/'run_r51_bridge_controlled.py',binary,formal,
        'bridge_'+identity+'_formal','--formal-gates',control/'formal_gates.json'])
    formal_state=json.loads((formal/'run_state.json').read_text())
    assert formal_state['returncode']==0 and formal_state['status']=='EXITED'
    assert formal_state.get('committed_ar_epochs')==120
    assert formal_state.get('last_committed_epoch')=='2024-07-17 02:59:30'
    state.update(status='COMPLETED_FORMAL_HOUR',formal_ar_started=True,completed_unix=time.time(),
        committed_ar_epochs=120,last_committed_epoch=formal_state['last_committed_epoch'],
        note='Completed calculation, not an independent absolute-product-accuracy certification')
except BaseException as exception:
    state.update(status='STOPPED_NO_AUTOMATIC_RETRY',error=repr(exception),finished_unix=time.time())
    if (formal/'run_state.json').exists():
        formal_state=json.loads((formal/'run_state.json').read_text())
        state['formal_ar_started']='solver_pid' in formal_state
    save()
    raise
save()
print(json.dumps(state,indent=2),flush=True)
