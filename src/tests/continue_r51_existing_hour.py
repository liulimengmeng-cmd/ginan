"""Continue the existing offline experiment; never start/restart a PEA.

The old supervisor has no live cancellation API. Suspend only its two verified
Python controllers so their 02:38 SIGINT and duplicate-hour launch cannot run.
Keep the actual solver, its timeout, input namespace and numerical policy intact.
The replacement observes the original launcher's exit status without reaping it.
No completion or scientific acceptance is inferred from a live process.
"""
from pathlib import Path
from datetime import datetime,timedelta
import hashlib
import json
import os
import shutil
import signal
import stat
import subprocess
import sys
import tempfile
import time

BASE=Path('/home/rx/GINAN/bridge_continuity_20261020')

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()

def process(pid):
    try:
        values=Path(f'/proc/{pid}/stat').read_text().rsplit(')',1)[1].split()
        return dict(pid=pid,state=values[0],ppid=int(values[1]),start_ticks=int(values[19]),
            exit_status=int(values[49]),cmdline=Path(f'/proc/{pid}/cmdline').read_bytes().replace(b'\0',b' ').decode())
    except FileNotFoundError:return None

def same_process(identity):
    current=process(identity['pid'])
    return current if current and current['start_ticks']==identity['start_ticks'] else None

def save(path,value):
    temporary=path.with_name(path.name+'.direct.tmp')
    temporary.write_text(json.dumps(value,indent=2))
    temporary.replace(path)

def self_test():
    # Exercise ownership, direct-file output, stop/resume, real zombie status.
    with tempfile.TemporaryDirectory(prefix='r51-supervisor-') as directory:
        with (Path(directory)/'output.log').open('wb') as output:
            child=subprocess.Popen([sys.executable,'-c','import time;time.sleep(60)'],stdout=output,stderr=output)
            try:
                identity=process(child.pid)
                assert identity and same_process(identity)
                assert stat.S_ISREG(os.stat(f'/proc/{child.pid}/fd/1').st_mode)
                os.kill(child.pid,signal.SIGSTOP)
                for unused in range(100):
                    if same_process(identity)['state']=='T':break
                    time.sleep(.01)
                assert same_process(identity)['state']=='T'
                os.kill(child.pid,signal.SIGCONT)
                child.terminate()
                for unused in range(100):
                    result=same_process(identity)
                    if result and result['state']=='Z':break
                    time.sleep(.01)
                assert result['state']=='Z' and os.waitstatus_to_exitcode(result['exit_status'])==-signal.SIGTERM
                assert not same_process({**identity,'start_ticks':identity['start_ticks']+1})
            finally:
                if child.poll() is None:os.kill(child.pid,signal.SIGCONT);child.kill()
                child.wait()
    print('PASS existing-hour supervisor ownership/direct-output/STOP/CONT/zombie-exit/PID-reuse tests')

def main():
    if sys.argv[1:]==['--self-test']:self_test();return
    assert len(sys.argv)==3 and os.geteuid()==0,'Usage: CONTROLLED_EVIDENCE BUILD (WSL root)'
    evidence=Path(sys.argv[1]).resolve(strict=True)
    build=Path(sys.argv[2]).resolve(strict=True)
    assert evidence==BASE/'controlled_verified_v2'
    workflow_path=BASE/'workflow_verified_v2.json'
    state_path=evidence/'run_state.json'
    old_workflow=json.loads(workflow_path.read_text())
    old_state=json.loads(state_path.read_text())
    assert old_workflow['status']=='CONTROLLED_EVENT_REPLAY' and old_state['status']=='RUNNING'
    assert old_state['stop_at']=='02:38:00' and not old_state['events']
    assert old_state['last_committed_epoch']<'2024-07-17 02:38:00'
    solver=process(old_state['solver_pid'])
    launcher=process(old_state['launcher_pid'])
    controllers=[process(old_workflow['pid']),process(old_workflow['child_pid'])]
    assert solver and launcher and all(controllers)
    assert 'run_r51_bridge_workflow.py' in controllers[0]['cmdline']
    assert 'run_r51_bridge_controlled.py' in controllers[1]['cmdline']
    assert controllers[1]['ppid']==controllers[0]['pid'] and launcher['ppid']==controllers[1]['pid']
    assert '--stop-at 02:38:00' in controllers[1]['cmdline']
    assert Path(f"/proc/{solver['pid']}/exe").resolve()==Path(old_state['binary'])
    assert digest(f"/proc/{solver['pid']}/exe")==old_state['binary_sha256']
    runner=Path(old_state['run_root'])/'runner.log'
    for pid in (solver['pid'],launcher['pid']):
        for descriptor in (1,2):
            assert stat.S_ISREG(os.stat(f'/proc/{pid}/fd/{descriptor}').st_mode)
            assert Path(f'/proc/{pid}/fd/{descriptor}').resolve()==runner
    timeout=process(solver['ppid'])
    assert timeout and 'timeout --signal=TERM --kill-after=30s 14d ' in timeout['cmdline']
    config=Path('/home/rx/GINAN/frozen-r51-satellite-fix-20260918/config/case.yaml')
    pinned=subprocess.check_output(['nsenter','-t',str(solver['pid']),'-m','cat',str(config)],text=True)
    assert 'end_epoch: 2024-07-17 02:59:30' in pinned
    adoption=BASE/'direct_hour_verified_v2'
    adoption.mkdir(exist_ok=False)
    shutil.copy2(__file__,adoption/'supervisor_used.py')
    auditor=Path(__file__).with_name('audit_r51_bridge_controlled.py')
    shutil.copy2(auditor,adoption/'auditor_used.py')
    audit_script=adoption/'auditor_used.py'
    save(adoption/'original_workflow.json',old_workflow)
    save(adoption/'original_run_state.json',old_state)
    record=dict(status='PREPARING',pid=os.getpid(),solver=solver,launcher=launcher,controllers=controllers,
        supervisor_sha256=digest(adoption/'supervisor_used.py'),auditor_sha256=digest(audit_script),
        authorization='User requested that this existing run continue directly for the full hour',
        no_solver_restart=True,numerical_policy_unchanged=True,duplicate_formal_run_disabled=True,
        original_stop_at='02:38:00',end_epoch='2024-07-17 02:59:30',scientific_acceptance=False)
    save(adoption/'adoption_state.json',record)
    stopped=[];active=False
    try:
        for identity in controllers:
            assert same_process(identity)
            os.kill(identity['pid'],signal.SIGSTOP);stopped.append(identity)
        for unused in range(100):
            if all(same_process(identity) and same_process(identity)['state']=='T' for identity in controllers):break
            time.sleep(.01)
        assert all(same_process(identity)['state']=='T' for identity in controllers)
        assert same_process(solver) and same_process(solver)['state'] not in ('T','Z')
        state=json.loads(state_path.read_text())
        assert not state['events'] and state['last_committed_epoch']<'2024-07-17 02:38:00'
        save(adoption/'suspended_run_state.json',state)
        state.update(stage='full_hour_in_place_pending_acceptance',original_stop_at='02:38:00',stop_at=None,
            end_epoch='2024-07-17 02:59:30',supervisor_pid=os.getpid(),scientific_acceptance=False,
            qualification='USER_AUTHORIZED_OFFLINE_HOUR_PENDING_EVENT_ACCEPTANCE',
            execution_plan='ONE_CONTINUOUS_HOUR_NO_SECOND_RUN')
        state['events'].append(dict(event='USER_AUTHORIZED_CONTINUE_EXISTING_RUN_WITHOUT_RESTART',
            unix=time.time(),old_controllers_suspended=[item['pid'] for item in controllers]))
        workflow={**old_workflow,'status':'FULL_HOUR_IN_PLACE_RUNNING_PENDING_ACCEPTANCE',
            'pid':os.getpid(),'superseded_workflow_pid':old_workflow['pid'],'child_pid':solver['pid'],
            'full_hour_in_place':True,'formal_ar_started':False,'duplicate_formal_run_disabled':True,
            'execution_plan':'ONE_CONTINUOUS_HOUR_NO_SECOND_RUN',
            'scientific_acceptance':False,'adoption_evidence':str(adoption)}
        save(state_path,state);save(workflow_path,workflow)
        record.update(status='ACTIVE',activated_unix=time.time())
        save(adoption/'adoption_state.json',record);active=True
        print(json.dumps(dict(status='ACTIVE',solver_pid=solver['pid'],unchanged_solver_start_ticks=solver['start_ticks'],
            end_epoch=state['end_epoch'],suspended_controllers=[item['pid'] for item in controllers])),flush=True)
        while True:
            current_launcher=same_process(launcher)
            assert current_launcher,'Original launcher disappeared without recorded exit status'
            if current_launcher['state']=='Z':break
            assert all(same_process(identity) and same_process(identity)['state']=='T' for identity in controllers), \
                'Old stop/restart controller was unexpectedly resumed'
            live=same_process(solver)
            commits=[]
            for path in (Path(state['run_root'])/'outputs/zhang_internal_products.csv.epochs').glob('*/COMMIT.json'):
                item=json.loads(path.read_text())
                if item.get('ar') is True:commits.append(item['epoch'])
            commits.sort()
            state.update(committed_ar_epochs=len(commits),last_committed_epoch=commits[-1] if commits else None,
                solver_still_running=bool(live and live['state']!='Z'),elapsed_seconds=time.time()-old_workflow['started_unix'])
            save(state_path,state)
            if time.time()-old_workflow['started_unix']>36*3600:
                raise TimeoutError('36-hour supervisory limit reached; solver not auto-killed or restarted')
            time.sleep(5)
        code=os.waitstatus_to_exitcode(current_launcher['exit_status'])
        live=same_process(solver)
        assert not live or live['state']=='Z','Launcher exited while the owned solver is still running'
        # Final rescan includes a commit emitted immediately before normal exit.
        commits=sorted(json.loads(p.read_text())['epoch'] for p in
            (Path(state['run_root'])/'outputs/zhang_internal_products.csv.epochs').glob('*/COMMIT.json')
            if json.loads(p.read_text()).get('ar') is True)
        expected=[(datetime(2024,7,17,2)+timedelta(seconds=30*i)).strftime('%Y-%m-%d %H:%M:%S') for i in range(120)]
        complete=code==0 and commits==expected
        state.update(status='EXITED' if code==0 else 'FAILED',returncode=code,solver_still_running=False,
            committed_ar_epochs=len(commits),last_committed_epoch=commits[-1] if commits else None,
            full_hour_calculation_complete=complete,scientific_acceptance=False)
        assert digest(state['checkpoint'])==state['checkpoint_sha256'],'Original checkpoint changed'
        assert digest(Path(state['checkpoint']).with_name('checkpoint_manifest.json'))==state['checkpoint_manifest_sha256']
        save(state_path,state)
        workflow.update(status='COMPLETED_FULL_HOUR_PENDING_ACCEPTANCE' if complete else 'FAILED_OR_INCOMPLETE_HOUR',
            completed_unix=time.time(),committed_ar_epochs=len(commits),full_hour_calculation_complete=complete)
        save(workflow_path,workflow)
        record.update(status=workflow['status'],launcher_exit_status=current_launcher['exit_status'],
            native_returncode=code,full_hour_calculation_complete=complete)
        save(adoption/'adoption_state.json',record)
        # Cleanup only verified obsolete Python controllers AFTER solver exit and durable evidence.
        for identity in controllers:
            if same_process(identity):os.kill(identity['pid'],signal.SIGKILL)
        if complete:
            with (adoption/'full_hour_audit.log').open('w') as output:
                result=subprocess.run([sys.executable,str(audit_script),str(evidence),str(build),'--completed-hour'],
                    stdout=output,stderr=subprocess.STDOUT)
            workflow.update(status='COMPLETED_FULL_HOUR_ACCEPTANCE_PASSED' if result.returncode==0 else
                'COMPLETED_FULL_HOUR_ACCEPTANCE_FAILED',acceptance_passed=result.returncode==0)
            # This acceptance scope is event/covariance/integer consistency, not independent absolute accuracy.
            save(workflow_path,workflow)
            state['event_consistency_acceptance_passed']=result.returncode==0
            save(state_path,state)
            record.update(status=workflow['status'],audit_returncode=result.returncode)
            save(adoption/'adoption_state.json',record)
    except BaseException as exception:
        record.update(status='SUPERVISOR_FAILED_SOLVER_NOT_RESTARTED',error=repr(exception))
        save(adoption/'adoption_state.json',record)
        if not active:
            for identity in reversed(stopped):
                if same_process(identity):os.kill(identity['pid'],signal.SIGCONT)
        # After takeover, never resume the 02:38 stopper: it could interrupt the run.
        raise

if __name__=='__main__':main()
