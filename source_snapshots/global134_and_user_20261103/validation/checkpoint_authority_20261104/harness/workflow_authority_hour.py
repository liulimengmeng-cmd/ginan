"""One authorised finite workflow. No retries, resource stops or server access."""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import secrets
import re
import signal
import subprocess
import sys
import time
from types import SimpleNamespace

HERE=Path(__file__).resolve().parent
ART=Path('/mnt/c/Users/rx/Documents/GINAN/artifacts/checkpoint_authority_fix_20261104_01')
ROOT=Path('/home/rx/GINAN/checkpoint_authority_fix_20261104_01')
DIRECTORY=ART/'workflow_01'
SCHEMA='CHECKPOINT_AUTHORITY_ONCE_WORKFLOW_PLAN_V1'
GATE_SCHEMA='CHECKPOINT_AUTHORITY_ONCE_WORKFLOW_GATE_V1'
RESULT_SCHEMA='CHECKPOINT_AUTHORITY_ONCE_WORKFLOW_RESULT_V1'
R=None

def require(ok,reason):
    if not ok:raise RuntimeError(reason)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1048576),b''):h.update(block)
    return h.hexdigest()

def identity(path):
    path=Path(path)
    return dict(path=str(path),bytes=path.stat().st_size,sha256=sha(path))

def load(path):
    def pairs(rows):
        out={}
        for key,value in rows:
            require(key not in out,'Duplicate JSON key');out[key]=value
        return out
    def number(value):
        actual=float(value);require(math.isfinite(actual),'Overflowed JSON number');return actual
    return json.loads(Path(path).read_text(),object_pairs_hook=pairs,parse_float=number,
        parse_constant=lambda v:(_ for _ in ()).throw(RuntimeError('Nonfinite JSON')))

def verify(pin):
    require(isinstance(pin,dict) and set(pin)=={'path','bytes','sha256'}
        and isinstance(pin.get('path'),str) and Path(pin['path']).is_absolute()
        and type(pin.get('bytes')) is int and pin['bytes']>0
        and isinstance(pin.get('sha256'),str) and re.fullmatch('[0-9a-f]{64}',pin['sha256'])
        and identity(pin['path'])==pin,
        'Pinned workflow dependency changed')

def write_new(path,value):
    with Path(path).open('x',encoding='utf-8',newline='\n') as stream:
        json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n');stream.flush();os.fsync(stream.fileno())

def atomic_status(value):
    path=DIRECTORY/'STATUS.json';temp=DIRECTORY/'STATUS.json.tmp'
    temp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');os.replace(temp,path)

def exclusive_failure(args,error,phase,**context):
    """Never grants canonical ownership or replaces an earlier failure file."""
    nonce=getattr(args,'nonce',None)
    tag=nonce if isinstance(nonce,str) and re.fullmatch('[0-9a-f]{32}',nonce) else 'unverified-'+hashlib.sha256(repr(nonce).encode()).hexdigest()[:24]
    actual=kernel();who=owner_only(actual) if actual is not None else None
    receipt=dict(schema='CHECKPOINT_AUTHORITY_WORKFLOW_EXCLUSIVE_FAILURE_V1',passed=False,phase=phase,
        nonce=nonce,owner=who,workflow_plan_sha256=getattr(args,'plan_sha256',None),error=repr(error),
        context=context,canonical_ownership_not_granted_by_this_receipt=True,
        unknown_processes_signalled=False,automatic_retry_performed=False,utc=datetime.now(timezone.utc).isoformat())
    path=DIRECTORY/(phase+'--'+tag+'--pid'+str(os.getpid())+'--birth'+str(who['birth_ticks'] if who else 'unavailable')+'.json')
    persisted=False
    try:
        require(DIRECTORY.is_dir() and not DIRECTORY.is_symlink(),'No prepared workflow evidence root')
        write_new(path,receipt);persisted=True
    except BaseException as secondary:
        receipt['persistence_error']=repr(secondary)
    print(json.dumps(dict(exclusive_workflow_failure=receipt,path=str(path),persisted=persisted),allow_nan=False),file=sys.stderr,flush=True)
    return receipt

def runner():
    global R
    if R is None:
        spec=importlib.util.spec_from_file_location('_authority_workflow_engine',HERE/'run_authority_hour.py')
        R=importlib.util.module_from_spec(spec);spec.loader.exec_module(R)
    return R

def nodes():
    return [dict(profile=p,action=a) for p in ('replay_e300','formal') for a in ('validation','formal','audit')]

def kernel(pid=None):
    pid=os.getpid() if pid is None else pid
    try:
        proc=Path('/proc')/str(pid);raw=(proc/'stat').read_text();fields=raw[raw.rfind(')')+2:].split()
        uid=int(next(line.split()[1] for line in (proc/'status').read_text().splitlines() if line.startswith('Uid:')))
        executable=os.readlink(proc/'exe');last=(proc/'stat').read_text();last=last[last.rfind(')')+2:].split()
        require(fields[19]==last[19],'PID changed during identity capture')
        return dict(pid=pid,birth_ticks=int(last[19]),sid=int(last[3]),uid=uid,
            executable=executable,state=last[0])
    except (FileNotFoundError,ProcessLookupError):return None

def same_owner(expected,actual):
    return isinstance(expected,dict) and isinstance(actual,dict) and all(
        type(expected.get(k)) is int and type(actual.get(k)) is int and expected[k]==actual[k]
        for k in ('pid','birth_ticks','sid','uid')) and bool(expected.get('executable')) and expected['executable']==actual.get('executable')

def owner_only(record):
    require(record is not None,'Missing process identity')
    return {key:record[key] for key in ('pid','birth_ticks','sid','uid','executable')}

def validate_contract(plan):
    require(plan.get('schema')==SCHEMA and plan.get('nodes')==nodes()
        and plan.get('original_source')=='/home/rx/GINAN/us20_global134_osb5h_20261031_v1/run/checkpoints/E29-20240717015930-e240-r32f626de9463'
        and plan.get('automatic_retry_authorised') is False and plan.get('resource_abort_thresholds_enabled') is False
        and plan.get('server_access_authorised') is False and plan.get('new_solver_started_by_prepare') is False
        and plan.get('formal_requires_actual_e300_restore') is True,'Changed finite conditional workflow')
    require(set(plan.get('execution_plans',{}))=={'replay_e300','formal'},'Both exact execution plans required')

def verify_plan(path,digest):
    require(Path(path)==DIRECTORY/'WORKFLOW_PLAN.json' and sha(path)==digest,'Exact workflow plan required')
    plan=load(path);validate_contract(plan)
    for pin in plan['pins']+list(plan['execution_plans'].values()):verify(pin)
    return plan

def validate_gate(gate,plan_sha):
    require(gate.get('schema')==GATE_SCHEMA and gate.get('approved') is True
        and gate.get('workflow_plan_sha256')==plan_sha and gate.get('nodes')==nodes()
        and gate.get('automatic_retry_authorised') is False
        and gate.get('allow_formal_only_after_all_replay_gates_pass') is True,
        'Exact root-reviewed conditional workflow gate required')

def completed_prefix(completed,index):
    require(len(completed)==index,'Completed workflow prefix has gaps')
    for position,row in enumerate(completed):
        require(row.get('node')==nodes()[position] and row.get('passed') is True,
            'Only actual passed ordered predecessors can authorise the next node')

def verify_completed_prefix(plan,completed,index):
    completed_prefix(completed,index)
    for row in completed:
        verify(row['result']);verify(row['launch'])
        result=load(row['result']['path']);launch=load(row['launch']['path']);node=row['node']
        expected_sha=plan['execution_plans'][node['profile']]['sha256']
        validate_result(result,node,expected_sha)
        require(launch.get('execution_plan_sha256')==expected_sha and launch.get('action')==node['action']
            and result.get('action_nonce')==launch.get('action_nonce')
            and same_owner(result.get('controller_identity'),launch.get('controller')),
            'Saved predecessor launch/result identity or nonce changed')
        proof=row.get('exact_child_waitpid',{});expected=launch['controller']
        require(proof.get('pid')==expected['pid'] and proof.get('birth_ticks')==expected['birth_ticks']
            and proof.get('exit_code')==0 and not same_owner(expected,kernel(expected['pid'])),
            'Earlier exact controller exit is not established before next action')

def validate_result(result,node,execution_sha):
    require(result.get('schema')=='CHECKPOINT_AUTHORITY_AR_ACTION_RESULT_V1'
        and result.get('passed') is True and result.get('profile')==node['profile']
        and result.get('action')==node['action'] and result.get('execution_plan_sha256')==execution_sha
        and result.get('original_source_objects_configuration_inputs_unchanged') is True,
        'Stage did not pass exact current-plan postflight')
    count=60 if node['profile']=='replay_e300' else 120
    if node['action']=='validation':
        require(result.get('actual_fresh_PEA_restores')==result.get('all_nine_section_audits')==2
            and len(result.get('stages',[]))==2 and result.get('finite_data_run_started') is False,
            'Both source whole-process restores and independent audits required')
    else:
        require(result.get('complete_epochs')==result.get('complete_FLOAT_epochs')==result.get('AR_enabled_epochs')==count,
            'All finite epoch boundaries required')
        if node['action']=='audit':
            captures=1 if node['profile']=='replay_e300' else 2
            require(result.get('atomic_product_commits')==count and len(result.get('products',[]))==count
                and len(result.get('regression_product_events',[]))==5
                and all(row.get('passed') is True for row in result['regression_product_events'])
                and result.get('all_captured_checkpoints_actual_roundtrip_tested') is True
                and len(result.get('actual_checkpoint_restores',[]))==captures
                and set(result.get('checkpoints',{}))==({'300'} if captures==1 else {'300','360'}),
                'All actual checkpoint restores and atomic product audits required')

def prepare(args):
    R=runner();plans={}
    require(sha(args.workflow_test_receipt)==args.workflow_test_sha256,'Exact workflow pure-test receipt required')
    pure=load(args.workflow_test_receipt)
    require(pure.get('passed') is True and pure.get('tests_run',0)>=14 and pure.get('solver_processes_started')==0
        and pure.get('workflow')==identity(Path(__file__).resolve())
        and pure.get('test_source')==identity(HERE/'test_authority_workflow.py'),
        'Current workflow/test source has not passed its actual pure tests')
    for profile,path,digest in (('replay_e300',args.replay_plan,args.replay_sha256),('formal',args.formal_plan,args.formal_sha256)):
        engine=R.configure(profile);engine.local_only();value=engine.verify_plan(path,digest,full_inputs=False)
        require(value['science']['profile']==profile,'Execution-plan profile mismatch')
        require(not any((ART/('controller_'+profile+'_01')/(action+'_REQUESTED.json')).exists()
            for action in ('validation','formal','audit')),'Workflow cannot adopt an already launched action')
        plans[profile]=(identity(path),value)
    require(plans['replay_e300'][1]['qualification']==plans['formal'][1]['qualification']
        and plans['replay_e300'][1]['science']['source']==plans['formal'][1]['science']['source']
        and plans['replay_e300'][1]['science']['target']['binary_sha256']==plans['formal'][1]['science']['target']['binary_sha256'],
        'Both profiles must use the same qualified new candidate and original source')
    require(not DIRECTORY.exists(),'Fresh single-use workflow directory required')
    DIRECTORY.mkdir()
    plan=dict(schema=SCHEMA,nodes=nodes(),execution_plans={key:value[0] for key,value in plans.items()},
        original_source=plans['formal'][1]['science']['source']['directory'],
        candidate=plans['formal'][1]['executables'][0],qualification=plans['formal'][1]['qualification'],
        automatic_retry_authorised=False,resource_abort_thresholds_enabled=False,server_access_authorised=False,
        new_solver_started_by_prepare=False,formal_requires_actual_e300_restore=True,
        pins=[identity(path) for path in (Path(__file__).resolve(),HERE/'run_authority_hour.py',
            HERE/'prepare_authority_hour.py',HERE/'test_authority_hour.py',HERE/'test_authority_workflow.py',args.workflow_test_receipt)])
    validate_contract(plan);write_new(DIRECTORY/'WORKFLOW_PLAN.json',plan)
    print(json.dumps(dict(prepared=True,plan=identity(DIRECTORY/'WORKFLOW_PLAN.json'),solver_started=False)))

def launch(args):
    plan=verify_plan(args.plan,args.plan_sha256);verify(dict(path=str(args.gate),bytes=Path(args.gate).stat().st_size,sha256=args.gate_sha256))
    validate_gate(load(args.gate),args.plan_sha256)
    runner().configure('replay_e300').local_only()
    require(not any((DIRECTORY/name).exists() for name in ('REQUESTED.json','LAUNCH.json','STARTED.json','RESULT.json')),
        'Single-use workflow; no restart or adoption')
    require(ROOT.is_dir() and not ROOT.is_symlink() and ROOT.resolve(strict=True)==ROOT,'Canonical local candidate root required')
    nonce=secrets.token_hex(16);launcher=owner_only(kernel())
    command=['/usr/bin/python3','-I','-S','-B',str(Path(__file__).resolve()),'--run','--plan',str(args.plan),
        '--plan-sha256',args.plan_sha256,'--gate',str(args.gate),'--gate-sha256',args.gate_sha256,'--nonce',nonce]
    write_new(DIRECTORY/'REQUESTED.json',dict(command=command,nonce=nonce,launcher=launcher,plan_sha256=args.plan_sha256,
        gate=identity(args.gate),automatic_retry_authorised=False))
    child=None;descriptor=None
    try:
        with (DIRECTORY/'controller.log').open('xb') as output:
            child=subprocess.Popen(command,cwd=ROOT,env=runner().M.native_environment({}),stdout=output,
                stderr=subprocess.STDOUT,start_new_session=True)
        descriptor=os.pidfd_open(child.pid,0);actual=kernel(child.pid)
        require(actual and actual['sid']==child.pid and actual['uid']==1000 and
            actual['executable']==str(Path('/usr/bin/python3').resolve()),'Workflow kernel identity not proven')
        owner=owner_only(actual)
        write_new(DIRECTORY/'LAUNCH.json',dict(command=command,nonce=nonce,controller=owner,
            launcher=launcher,plan_sha256=args.plan_sha256,gate=identity(args.gate)))
        print(json.dumps(dict(launched=True,controller=owner,launch=identity(DIRECTORY/'LAUNCH.json'))))
    except BaseException:
        # No solver starts until this exact handshake exists. Signal only the
        # launcher-created process whose kernel handle was obtained above.
        if descriptor is not None:signal.pidfd_send_signal(descriptor,signal.SIGTERM,None,0)
        raise
    finally:
        if descriptor is not None:os.close(descriptor)

def wait_own_launch(args,owner):
    require(isinstance(args.nonce,str) and len(args.nonce)==32 and all(c in '0123456789abcdef' for c in args.nonce),
        'Exact one-use workflow nonce required')
    deadline=time.monotonic()+30
    while not (DIRECTORY/'LAUNCH.json').exists():
        require(time.monotonic()<deadline,'Workflow launch identity handshake did not complete');time.sleep(.05)
    receipt=load(DIRECTORY/'LAUNCH.json')
    require(receipt.get('nonce')==args.nonce and receipt.get('plan_sha256')==args.plan_sha256
        and same_owner(receipt.get('controller'),owner),'Workflow cannot claim another launch')
    require(receipt['gate']==identity(args.gate) and sha(args.gate)==args.gate_sha256,'Launch gate changed')

def child_status(profile,action):
    path=ROOT/('controller_'+profile+'_01')/action/'status.json'
    if not path.exists():return None
    value=load(path)
    return {key:value.get(key) for key in ('utc','stage','actual_phase','active_epoch','last_complete_ar_epoch',
        'complete_AR_epochs','errors')}

def action_gate(plan,index,args,owner,completed):
    verify_completed_prefix(plan,completed,index);node=nodes()[index];profile,action=node['profile'],node['action']
    execution=plan['execution_plans'][profile];evidence=ART/('controller_'+profile+'_01')
    value=dict(schema=runner().GATE_SCHEMA,action=action,profile=profile,approved=True,
        execution_plan_sha256=execution['sha256'],automatic_retry_authorised=False,
        conditional_workflow_plan=identity(args.plan),conditional_workflow_gate=identity(args.gate),
        conditional_workflow_nonce=args.nonce,conditional_workflow_owner=owner,node_index=index,
        passed_ordered_predecessors=[row['result'] for row in completed])
    if action in ('formal','audit'):
        previous=evidence/(('validation' if action=='formal' else 'formal')+'_RESULT.json')
        value['predecessor_result']=identity(previous)
    if profile=='formal':value['replay_acceptance']=identity(ART/'controller_replay_e300_01/audit_RESULT.json')
    target=DIRECTORY/(f'{index+1:02d}_{profile}_{action}_GATE.json');write_new(target,value)
    return target,identity(target)

def wait_child(plan,index,owner,args,completed,launch):
    node=nodes()[index];pid=launch['controller']['pid'];expected=launch['controller']
    execution=plan['execution_plans'][node['profile']]
    evidence=ART/('controller_'+node['profile']+'_01');result_path=evidence/(node['action']+'_RESULT.json')
    while True:
        actual=kernel(pid)
        if actual is not None:require(same_owner(expected,actual),'Owned controller changed identity before exit')
        # This controller is a direct child of this workflow. Reap the exact
        # child before any next Popen; an absent result is not success.
        waited,status=os.waitpid(pid,os.WNOHANG)
        if waited==pid:
            exit_code=os.waitstatus_to_exitcode(status)
            require(exit_code==0,'Owned stage controller exited with '+str(exit_code))
            require(result_path.exists(),'Stage controller exited without its actual result')
            result=load(result_path)
            require(result.get('action_nonce')==launch['action_nonce'] and
                same_owner(result.get('controller_identity'),expected),'Stage result owner/nonce mismatch')
            validate_result(result,node,execution['sha256'])
            return dict(node=node,passed=True,result=identity(result_path),launch=identity(evidence/(node['action']+'_LAUNCH.json')),
                exact_child_waitpid=dict(pid=pid,birth_ticks=expected['birth_ticks'],exit_code=exit_code),
                original_controller_absent_after_wait=not same_owner(expected,kernel(pid)))
        require(waited==0 and same_owner(expected,actual),'Active stage controller kernel identity changed')
        atomic_status(dict(state='RUNNING',owner=owner,workflow_plan_sha256=args.plan_sha256,
            node_index=index,node=node,completed_nodes=len(completed),active_controller=expected,
            native_status=child_status(node['profile'],node['action']),utc=datetime.now(timezone.utc).isoformat()))
        time.sleep(2)

def run(args):
    owner=owner_only(kernel());wait_own_launch(args,owner)
    require(not (DIRECTORY/'STARTED.json').exists() and not (DIRECTORY/'RESULT.json').exists(),
        'Workflow cannot resume an earlier attempt')
    write_new(DIRECTORY/'STARTED.json',dict(owner=owner,nonce=args.nonce,plan_sha256=args.plan_sha256,
        utc=datetime.now(timezone.utc).isoformat()))
    completed=[];index=-1;active=None
    outcome=dict(schema=RESULT_SCHEMA,passed=False,workflow_plan_sha256=args.plan_sha256,owner=owner,
        nonce=args.nonce,completed=completed,automatic_retry_performed=False,scientific_accuracy_accepted=False)
    try:
        plan=verify_plan(args.plan,args.plan_sha256);validate_gate(load(args.gate),args.plan_sha256)
        require(sha(args.gate)==args.gate_sha256,'Workflow gate SHA changed')
        for index,node in enumerate(nodes()):
            verify_plan(args.plan,args.plan_sha256);verify_completed_prefix(plan,completed,index)
            engine=runner().configure(node['profile']);engine.local_only()
            gate,pin=action_gate(plan,index,args,owner,completed)
            execution=plan['execution_plans'][node['profile']]
            engine.launch(SimpleNamespace(launch=node['action'],execution_plan=Path(execution['path']),
                plan_sha256=execution['sha256'],gate=gate,gate_sha256=pin['sha256']))
            active=load(ART/('controller_'+node['profile']+'_01')/(node['action']+'_LAUNCH.json'))
            completed.append(wait_child(plan,index,owner,args,completed,active));active=None
        final=ART/'controller_formal_01/hour_result.json'
        final_result=load(final)
        require(final_result.get('passed') is True and final_result.get('hour_complete') is True
            and final_result.get('complete_epochs')==final_result.get('atomic_product_commits')==120,
            'Formal hour independent final acceptance is incomplete')
        outcome.update(passed=True,formal_hour_result=identity(final),complete_new_AR_epochs=120,
            e300_replay_restored_and_audited=True,formal_e300_e360_restored_and_audited=True)
    except BaseException as error:
        outcome.update(error=repr(error),failed_node_index=index,
            failed_node=nodes()[index] if index>=0 else None,active_controller=active,
            no_further_actions_started=True,
            active_owned_stage_may_continue_if_workflow_infrastructure_failed=active is not None)
        raise
    finally:
        outcome['completed']=completed;outcome['utc']=datetime.now(timezone.utc).isoformat()
        try:
            write_new(DIRECTORY/'RESULT.json',outcome)
            atomic_status(dict(state='PASSED' if outcome['passed'] else 'FAILED',result=identity(DIRECTORY/'RESULT.json'),
                completed_nodes=len(completed),utc=outcome['utc']))
        except BaseException as finalization_error:
            exclusive_failure(args,finalization_error,'FINALIZATION',original_action_error=outcome.get('error'),
                completed_nodes=len(completed))
            if 'error' not in outcome:raise

def main():
    p=argparse.ArgumentParser(description=__doc__);mode=p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare',action='store_true');mode.add_argument('--launch',action='store_true');mode.add_argument('--run',action='store_true')
    for name in ('replay-plan','formal-plan','plan','gate','workflow-test-receipt'):p.add_argument('--'+name,type=Path)
    for name in ('replay-sha256','formal-sha256','plan-sha256','gate-sha256','nonce','workflow-test-sha256'):p.add_argument('--'+name)
    args=p.parse_args()
    if args.prepare:prepare(args)
    elif args.launch:launch(args)
    else:
        try:run(args)
        except BaseException as error:
            if not (DIRECTORY/'RESULT.json').exists():
                exclusive_failure(args,error,'WITHOUT_CANONICAL_RESULT',started_exists=(DIRECTORY/'STARTED.json').exists())
            raise

if __name__=='__main__':main()
