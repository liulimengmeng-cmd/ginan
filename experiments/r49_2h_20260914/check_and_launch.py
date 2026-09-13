"""One bounded heartbeat check. Launch once only after a verified parent exit."""
from pathlib import Path
import argparse, datetime, fcntl, hashlib, json, re, shutil, subprocess

W = Path(__file__).resolve().parent

def read_json(path):
    try: return json.loads(path.read_text())
    except (FileNotFoundError, json.JSONDecodeError): return None

def completion(audit, spec):
    state = audit / (spec['case']+'.state.txt')
    text = state.read_text() if state.exists() else ''
    exits = re.findall(r'RUNNER_EXIT status=(\d+)',text)
    if not exits: return 'PENDING'
    if exits[-1]!='0': return 'FAILED'
    result = read_json(audit/'r49_final_verification.json')
    if not result: return 'FAILED'
    if (result.get('case') != spec['case'] or result.get('pea_status') != 0
        or result.get('runtime_pass') is not True or result.get('failures') or result.get('alarms')
        or len(result.get('epochs',[])) != spec['observations']['expected_epochs']):
        return 'FAILED'
    start=datetime.datetime.fromisoformat(spec['observations']['start']).replace(tzinfo=datetime.timezone.utc)
    end=datetime.datetime.fromisoformat(spec['observations']['end']).replace(tzinfo=datetime.timezone.utc)
    expected=list(range(int(start.timestamp()),int(end.timestamp())+1,spec['observations']['interval_s']))
    if [r['epoch'] for r in result['epochs']] != expected: return 'FAILED'
    return 'PASSED'

def processes():
    rows=[]
    for p in Path('/proc').iterdir():
        if not p.name.isdigit():continue
        try:
            args=(p/'cmdline').read_bytes().split(b'\0')
            rows.append((int(p.name), [a.decode(errors='replace') for a in args if a]))
        except (FileNotFoundError, PermissionError, ProcessLookupError):pass
    return rows

def check(launch=False):
    spec=read_json(W/'experiment.json')
    parent=Path(spec['continuation']['parent_audit'])
    pspec=read_json(parent/'experiment.json')
    rows=processes()
    def active(case):return [pid for pid,args in rows if case in args]
    if (W/'launch_claim.json').exists():
        outcome=completion(W,spec)
        if outcome!='PENDING':return {'status':'NEXT_'+outcome}
        running=active(spec['case'])
        runners=[pid for pid,args in rows if str(W/'run_r49_2h.sh') in args]
        return {'status':'NEXT_RUNNING' if running or runners else 'NEXT_NEEDS_ATTENTION', 'solver_pids':running,'runner_pids':runners}
    parent_status=completion(parent,pspec)
    if parent_status=='FAILED':return {'status':'PARENT_FAILED','audit':str(parent)}
    if parent_status=='PENDING':
        running=active(pspec['case'])
        runners=[pid for pid,args in rows if str(parent/'run_r49_30m.sh') in args]
        return {'status':'PARENT_RUNNING' if running or runners else 'PARENT_NEEDS_ATTENTION','solver_pids':running,'runner_pids':runners}
    pea=[pid for pid,args in rows if args and Path(args[0]).name=='pea']
    if pea:return {'status':'WAIT_OTHER_SOLVER','pids':pea}
    if not launch:return {'status':'READY_CHECK_ONLY'}
    # Never retry an ambiguous launch. Exclusive claim also survives a host restart.
    for line in (W/'PREPARED_SHA256SUMS').read_text().splitlines():
        digest,name=line.split('  ',1)
        if hashlib.sha256((W/name).read_bytes()).hexdigest()!=digest:
            return {'status':'PREPARED_HASH_FAILURE','file':name}
    out=Path('/mnt/d/GINAN_R20/inputData/outputs')/spec['case']
    for path in [out,Path(spec['checkpoints']['directory']),Path(spec['r49_policy']['root_snapshot_directory']),W/'launch.log',W/(spec['case']+'.runner.log'),W/(spec['case']+'.state.txt')]:
        if path.exists():return {'status':'OUTPUT_COLLISION','path':str(path)}
    free=shutil.disk_usage('/mnt/d/GINAN_R20').free
    if free<100*1024**3:return {'status':'INSUFFICIENT_DISK','free_bytes':free}
    receipt={'case':spec['case'],'claimed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'parent_verified':True,'expected_epochs':241}
    with (W/'launch_claim.json').open('x') as f:json.dump(receipt,f,indent=2)
    with (W/'launch.log').open('xb') as log:
        proc=subprocess.Popen(['bash',str(W/'run_r49_2h.sh')],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,close_fds=True)
    receipt['runner_pid']=proc.pid
    (W/'launch_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return {'status':'NEXT_LAUNCHED','receipt':receipt,'note':'Runner dispatched; verify solver startup separately.'}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--launch-if-ready',action='store_true');args=parser.parse_args()
    with (W/'heartbeat_check.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        result=check(args.launch_if_ready)
        result['checked_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
        (W/'heartbeat_latest.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result,indent=2))
