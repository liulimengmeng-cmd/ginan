"""Six user-approved 13:10--14:10 GPST cases, using the frozen solver/model.

No provider rerun, candidate preregistration, or source modification. Each case
starts in a fresh process and preserves all 120 requested posterior audit epochs.
"""
from pathlib import Path
import concurrent.futures as cf
import copy,datetime as dt,fcntl,io,json,math,os,shlex,shutil,struct,subprocess,sys,tarfile,time,traceback
import numpy as np
import yaml

R=Path('/home/rx/GINAN/us20_reset_gec_20261029_v1')
C=R/'candidates/preregistered_fork_014';D=C/'daily_direct_20261029_v1'
B=Path('/home/rx/GINAN/us20_osb_event_1327_1h_20261029_v1')
H=Path('/home/rx/GINAN/us20_hourly_restart_ar_20261029_v1')
PY='/home/rx/GINAN/france_users_ppp_gapfix_20261026_v1/plot_env/bin/python'
sys.path.insert(0,str(C))
import audit_daily_results_parallel_v3 as a
import run_direct_daily_7p5_v1 as base
sys.path.insert(0,str(base.U/'scripts'))
import calendar_accounting_v1 as accounting
import run_users as parser
START=dt.datetime(2026,10,1,13,10);END=dt.datetime(2026,10,1,14,9,30)
EVENT=dt.datetime(2026,10,1,13,27);N=120
load=base.load;save=base.save;sha=base.sha

def state(stage,**kw):save(B/'workflow_state.json',dict(state=stage,controller_pid=os.getpid(),updated_utc=base.now(),**kw))
def epochs():return [str(START+dt.timedelta(seconds=30*i)) for i in range(N)]
def matrix(p):
    data=p.read_bytes();n,m=struct.unpack('<QQ',data[:16]);assert len(data)==16+8*n*m
    return np.frombuffer(data[16:],dtype='<f8').reshape((n,m),order='F').copy()
def keys(p):
    result=[]
    for i,line in enumerate(p.read_text().splitlines()):
        f=shlex.split(line);assert int(f[0])==i
        result.append((int(f[1]),f[2],f[3],int(f[4]),f[5]))
    return result

def prepare(approval_path):
    assert not B.exists(),'Inspect existing root/processes; do not duplicate this window'
    approved=load(approval_path)
    assert approved['approved_window']==[str(START),'2026-10-01 14:10:00'] and approved['mode']=='AR'
    existing=load(H/'experiment.json')
    assert not any(x['start_epoch']==str(START) and x['end_epoch']==str(END) for x in existing['cases'])
    assert len(os.sched_getaffinity(0))>=12 and shutil.disk_usage(B.parent).free>140*1024**3
    manifest=load(D/'direct_run_manifest.json');assert manifest['ar_elevation_mask_deg']==7.5
    assert sha(manifest['binary'])==manifest['binary_sha256']=='fe75641d62a7abb9477f0c7d2b541e491ebc209b5f6dad7817ebdb454d90f9a8'
    products=R/'products/server_full_day_2880';export=load(products/'export.json')
    assert export['passed'] and export['epochs']==2880
    assert load(R/'cases/server_full_day_2880/full_day_acceptance.json')['passed']
    for name,v in export['files'].items():assert sha(products/name)==v['sha256']
    # Only prepare this already-authorized window; no preliminary solver jobs.
    B.mkdir();(B/'scripts').mkdir();shutil.copy2(__file__,B/'scripts'/Path(__file__).name)
    shutil.copy2(approval_path,B/'user_approval.json')
    reference=base.U/'reference_evaluation_only.json';shutil.copy2(reference,B/reference.name)
    rows=[];observations={}
    for item in [x for x in manifest['cases'] if x['mode']=='AR']:
        assert sha(item['config'])==item['config_sha256'] and sha(item['environment'])==item['environment_sha256']
        observations.setdefault(item['observation'],sha(item['observation']))
        assert observations[item['observation']]==item['observation_sha256_from_frozen_experiment']
        original=yaml.safe_load(Path(item['config']).read_text());oldenv=load(item['environment'])
        assert original['processing_options']['ambiguity_resolution']['elevation_mask']==7.5
        assert oldenv['ZHANG_USER_AR_SEARCH']=='GONSIS_UDUC' and oldenv['ZHANG_USER_AR_FEEDBACK']=='EPOCH_ONLY'
        assert oldenv['ZHANG_US20_EVENT_PROTOCOL']==oldenv['ZHANG_US20_RESET_CHECKPOINT']=='0'
        out=B/'cases'/item['policy']/(item['station']+'_AR');out.mkdir(parents=True)
        config=copy.deepcopy(original);config['outputs']['outputs_root']=str(out)
        config['outputs']['metadata']['config_description']='US20_OSB1327_1H_'+item['station']+'_'+item['policy']
        config['processing_options']['epoch_control'].update(start_epoch=str(START),end_epoch=str(END),max_epochs=N)
        opt=config['processing_options']['gnss_general']['zhang_pppar']
        assert opt['checkpoint_restore_path']=='' and not opt['deterministic_checkpoint']
        opt.update(checkpoint_output_directory=str(out/'checkpoints'),checkpoint_runtime_id='OSB1327_'+item['station']+'_'+item['policy'])
        env=copy.deepcopy(oldenv);env['ZHANG_RESET_AUDIT_DIRECTORY']=str(out/'state_audit')
        env['ZHANG_RESET_AUDIT_EPOCHS']=','.join(str(int((START-dt.datetime(1980,1,6)).total_seconds())+30*i) for i in range(N))
        env['ZHANG_USER_TROP_CONSUMPTION_LEDGER']=str(out/'trop_consumption_ledger')
        restored=copy.deepcopy(config);restored['outputs']=copy.deepcopy(original['outputs'])
        restored['processing_options']['epoch_control']=copy.deepcopy(original['processing_options']['epoch_control'])
        for k in ['checkpoint_output_directory','checkpoint_runtime_id']:restored['processing_options']['gnss_general']['zhang_pppar'][k]=original['processing_options']['gnss_general']['zhang_pppar'][k]
        assert restored==original
        restored_env=copy.deepcopy(env)
        for k in ['ZHANG_RESET_AUDIT_DIRECTORY','ZHANG_RESET_AUDIT_EPOCHS','ZHANG_USER_TROP_CONSUMPTION_LEDGER']:restored_env[k]=oldenv[k]
        assert restored_env==oldenv
        (out/'execution.yaml').write_text(yaml.safe_dump(config,sort_keys=False));save(out/'environment.json',env)
        row=dict(case=item['station']+'_AR',station=item['station'],policy=item['policy'],mode='AR',
            start_epoch=str(START),end_epoch=str(END),calendar_epochs=N,working_directory=str(out),
            config=str(out/'execution.yaml'),config_sha256=sha(out/'execution.yaml'),environment=str(out/'environment.json'),environment_sha256=sha(out/'environment.json'),
            observation=item['observation'],observation_sha256=observations[item['observation']],
            command=[manifest['binary'],'--rnx_inputs',item['observation'],'-y',str(out/'execution.yaml')],
            original_config=item['config'],original_config_sha256=item['config_sha256'],original_environment_sha256=item['environment_sha256'],
            scientific_settings_equal_except_window_and_output=True,posterior_audit_every_requested_epoch=True)
        rows.append(row);save(out/'runtime.json',dict(state='QUEUED',case=row['case'],station=row['station'],policy=row['policy']))
    assert len(rows)==6
    save(B/'experiment.json',dict(schema='US20_OSB1327_ONE_HOUR_SIX_AR_CASES_V1',created_utc=base.now(),cases=rows,
        start_GPST=str(START),stop_exclusive_GPST='2026-10-01 14:10:00',last_epoch_GPST=str(END),core_event_GPST=str(EVENT),
        binary=manifest['binary'],binary_sha256=manifest['binary_sha256'],reference_sha256=sha(reference),
        provider_export_sha256=sha(products/'export.json'),provider_files=export['files'],observation_hashes=observations,
        source_direct_manifest_sha256=sha(D/'direct_run_manifest.json'),ar_and_observation_elevation_deg=7.5,
        receiver_clock_handling=original['processing_options']['ppp_filter']['periodic_reset'],fresh_process_without_restore=True,
        AR_feedback='EPOCH_ONLY',six_cases_only=True,approval_sha256=sha(B/'user_approval.json'),
        frozen_statistics_sha256=a.FROZEN_ANALYSIS_SHA,
        limitations=['Short diagnostic window, not full-day validation','STATIC; frozen SOPAC parameter-model coordinates',
            'T/R are cumulative policies within this window, not a same-parent single-event intervention',
            'Partial antenna-frequency calibration','WHITE cross-epoch product error approximation','Fix fraction is not correct-integer probability']))
    state('PREPARED',cases=6,source_model_changed=False,new_scientific_solver_processes=0)
    print(json.dumps(dict(prepared=True,root=str(B),cases=6)),flush=True)

def postprocess(item):
    out=Path(item['working_directory']);rec=load(out/'runtime.json')
    assert load(out/'run.log.process.json')['returncode']==0
    accepted=accounting.validate_case_calendar(out,item,START,N)
    assert accepted['calendar_epochs']==N and accepted['native_PPP_status_epochs']==N and accepted['no_observation_epochs']==0
    reference=load(B/'reference_evaluation_only.json');av={x['epoch']:x for x in load(R/'products/server_full_day_2880/availability.json')}
    pts,sources=a.collect_points(out,item,reference['stations'][item['station']]['xyz_m'],av,epochs())
    policy_events=[];trop=[]
    for f in out.glob('Network-*.TRACE'):
        for line in f.open():
            if line.startswith('ZHANG_USER_DATUM_POLICY_EVENT '):policy_events.append(parser.event(line))
            elif line.startswith('ZHANG_USER_REGIONAL_TROP '):trop.append(parser.event(line))
    if item['policy']=='TRANSPORT_KNOWN':assert all(x['released_coordinates']=='0' for x in policy_events)
    save(out/'reset_event_receipts.json',dict(events=policy_events,trop=trop))
    ref=np.array(reference['stations'][item['station']]['xyz_m']);rot=np.array(a.M.rotation(ref.tolist()))
    core=[];samples=[]
    sample_epochs={str(EVENT-dt.timedelta(seconds=30)),str(EVENT),str(EVENT+dt.timedelta(seconds=30)),str(EVENT+dt.timedelta(minutes=10)),str(END)}
    for p in pts:
        if not p['fresh']:core.append(dict(p,enu_m=None,has_position=False));continue
        native=int((dt.datetime.fromisoformat(p['epoch'])-dt.datetime(1980,1,6)).total_seconds())
        q=out/'state_audit'/str(native)/'committed_posterior'
        kk=keys(q/'state_keys.txt');x=matrix(q/'x.f64').ravel();ix=[kk.index((3,item['station'],'',i,'')) for i in range(3)]
        enu=rot@(x[ix]-ref)
        core.append(dict(p,enu_m=enu.tolist(),xyz_m=x[ix].tolist(),horizontal_m=float(np.hypot(*enu[:2])),
            fixed=False,fixed_rank=0,fixed_rank_by_system={s:0 for s in a.SYSTEMS}))
        if p['epoch'] in sample_epochs:
            P=matrix(q/'P.f64');sigma=float(np.sqrt(rot[2]@P[np.ix_(ix,ix)]@rot[2]))
            ztd=next((i for i,k in enumerate(kk) if k[0]==21),None)
            samples.append(dict(epoch=p['epoch'],FLOAT_ENU_m=enu.tolist(),AR_output_ENU_m=p['enu_m'],sigma_FLOAT_U_m=sigma,
                ZTD_m=float(x[ztd]) if ztd is not None else None,
                rho_U_ZTD=float(rot[2]@P[ix,ztd]/(sigma*np.sqrt(P[ztd,ztd]))) if ztd is not None else None,
                source_hashes={n:sha(q/n) for n in ['state_keys.txt','x.f64','P.f64']}))
    stats=a.stats(pts,7200);restart=int((START-a.START).total_seconds());first=stats['convergence_primary']['first']
    stats['convergence_first_seconds_from_restart']=first['first_seconds_from_00']-restart if first else None
    stats['convergence_confirmation_seconds_from_restart']=first['confirmation_seconds_from_00']-restart if first else None
    stats['first_fix_seconds_from_restart']=stats['first_fix']['seconds']-restart if stats['first_fix'] else None
    save(out/'per_epoch.json',pts);save(out/'FLOAT_core_per_epoch.json',core)
    save(out/'native_FLOAT_AR_samples.json',samples)
    save(out/'statistics.json',dict(statistics=stats,FLOAT_core_accuracy=a.M.metrics(core),source_hashes=sources,
        event_receipts_at_1327=[e for e in policy_events if e['time']==str(EVENT)],native_posterior_samples=samples))
    rec.update(state='COMPLETE',finished_utc=base.now(),last_processed_calendar_epoch=N,last_processed_GPST=str(END),
        fresh_updates=accepted['fresh_updates'],acceptance_sha256=sha(out/'runtime_acceptance.json'))
    save(out/'runtime.json',rec);print('COMPLETE',item['station'],item['policy'],flush=True)

def one(item):
    out=Path(item['working_directory']);rec=load(out/'runtime.json')
    try:
        assert rec['state']=='QUEUED' and not (out/'run.log').exists()
        assert sha(item['config'])==item['config_sha256'] and sha(item['environment'])==item['environment_sha256']
        env={k:v for k,v in os.environ.items() if not k.startswith('ZHANG_')};env.update(load(item['environment']))
        with (out/'run.log').open('x') as log:
            begin=time.monotonic();p=subprocess.Popen(item['command'],cwd=out,env=env,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT)
            rec.update(state='RUNNING',solver_pid=p.pid,controller_pid=os.getpid(),started_utc=base.now());save(out/'runtime.json',rec)
            process=dict(state='RUNNING',solver_pid=p.pid,command=item['command'],started_utc=base.now());save(out/'run.log.process.json',process)
            print('STARTED',item['station'],item['policy'],'pid',p.pid,flush=True)
            while p.poll() is None:
                try:p.wait(timeout=20)
                except subprocess.TimeoutExpired:rec.update(base.progress(out/'run.log'),updated_utc=base.now());save(out/'runtime.json',rec)
            process.update(state='PROCESS_COMPLETE' if p.returncode==0 else 'FAILED',returncode=p.returncode,elapsed_seconds=time.monotonic()-begin,finished_utc=base.now());save(out/'run.log.process.json',process)
        assert p.returncode==0,('solver_failed',p.returncode)
        rec.update(state='POSTPROCESSING',solver_finished=True);save(out/'runtime.json',rec)
        postprocess(item)
    except BaseException:
        process=load(out/'run.log.process.json') if (out/'run.log.process.json').exists() else {}
        rec.update(state='POSTPROCESSING_FAILED' if process.get('returncode')==0 else 'FAILED',traceback=traceback.format_exc(),updated_utc=base.now())
        save(out/'runtime.json',rec);print('FAILED',item['station'],item['policy'],rec['traceback'],flush=True)

def summarize(exp):
    cases=[];pairs=[]
    points={}
    for item in exp['cases']:
        out=Path(item['working_directory']);assert load(out/'runtime.json')['state']=='COMPLETE'
        pts=load(out/'per_epoch.json');points[item['station'],item['policy']]=pts
        cases.append(dict(station=item['station'],policy=item['policy'],**load(out/'statistics.json'),
            before_event=a.stats([p for p in pts if p['epoch']<str(EVENT)],7200),after_event=a.stats([p for p in pts if p['epoch']>=str(EVENT)],7200)))
    for sta in ['MTA1','PKRD','P800']:
        left,right=[points[sta,p] for p in base.POLICIES]
        pairs.append(dict(station=sta,whole_hour=a.paired_stats(left,right,7200),
            after_event=a.paired_stats([p for p in left if p['epoch']>=str(EVENT)],[p for p in right if p['epoch']>=str(EVENT)],7200)))
    products=R/'products/server_full_day_2880'
    assert sha(exp['binary'])==exp['binary_sha256'] and sha(B/'reference_evaluation_only.json')==exp['reference_sha256']
    for name,v in exp['provider_files'].items():assert sha(products/name)==v['sha256']
    for item in exp['cases']:
        assert sha(item['config'])==item['config_sha256'] and sha(item['environment'])==item['environment_sha256']
    (B/'analysis').mkdir(exist_ok=True)
    save(B/'analysis/summary.json',dict(schema='OSB1327_ONE_HOUR_T_R_AR_RESULTS_V1',completed=True,created_utc=base.now(),
        experiment_sha256=sha(B/'experiment.json'),cases=cases,paired_T_R=pairs,limitations=exp['limitations']))
    files={'experiment.json':B/'experiment.json','reference_evaluation_only.json':B/'reference_evaluation_only.json',
        'user_approval.json':B/'user_approval.json','analysis/summary.json':B/'analysis/summary.json',
        'scripts/run_osb1327_one_hour_v1.py':B/'scripts/run_osb1327_one_hour_v1.py'}
    for item in exp['cases']:
        out=Path(item['working_directory'])
        for name in ['runtime.json','runtime_acceptance.json','run.log.process.json','run.log','calendar_status.json','event_receipts.json',
            'reset_event_receipts.json','per_epoch.json','FLOAT_core_per_epoch.json','native_FLOAT_AR_samples.json','statistics.json','execution.yaml','environment.json']:
            p=out/name;files[str(p.relative_to(B))]=p
        for p in out.glob('*.POS'):files[str(p.relative_to(B))]=p
    archive=B/'delivery_results.tar.gz';payload={}
    with tarfile.open(archive,'w:gz',compresslevel=3) as tar:
        for name,p in files.items():
            b=p.read_bytes();payload[name]=dict(bytes=len(b),sha256=sha(p));info=tarfile.TarInfo(name);info.size=len(b);tar.addfile(info,io.BytesIO(b))
        b=json.dumps(dict(files=payload),indent=2).encode();info=tarfile.TarInfo('payload_manifest.json');info.size=len(b);tar.addfile(info,io.BytesIO(b))
    state('COMPLETE_REMOTE_ANALYSIS_AWAITING_LOCAL_DELIVERY',completed_cases=6,summary=str(B/'analysis/summary.json'),
        summary_sha256=sha(B/'analysis/summary.json'),archive=str(archive),archive_sha256=sha(archive),archive_bytes=archive.stat().st_size)

def run():
    with (B/'controller.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        assert load(B/'workflow_state.json')['state']=='PREPARED'
        exp=load(B/'experiment.json');assert sha(exp['binary'])==exp['binary_sha256']
        with cf.ThreadPoolExecutor(max_workers=6) as pool:
            futures=[pool.submit(one,item) for item in exp['cases']]
            while not all(f.done() for f in futures):
                recs=[load(Path(x['working_directory'])/'runtime.json') for x in exp['cases']]
                state('RUNNING_OSB1327_ONE_HOUR',cases=6,completed=sum(r['state']=='COMPLETE' for r in recs),
                    running=sum(r['state']=='RUNNING' for r in recs),failed=sum('FAILED' in r['state'] for r in recs),max_concurrent_solvers=6)
                cf.wait(futures,timeout=20)
            for f in futures:f.result()
        recs=[load(Path(x['working_directory'])/'runtime.json') for x in exp['cases']]
        if any(r['state']!='COMPLETE' for r in recs):state('NEEDS_INSPECTION_NO_AUTOMATIC_RERUN',cases=recs);return
        summarize(exp)

if __name__=='__main__':
    if sys.argv[1:2]==['--prepare']:prepare(sys.argv[2])
    elif sys.argv[1:]==['--launch']:
        assert load(B/'workflow_state.json')['state']=='PREPARED'
        with (B/'controller.log').open('x') as log:
            p=subprocess.Popen([PY,str(B/'scripts/run_osb1327_one_hour_v1.py'),'--run'],cwd=B,
                env=dict(os.environ,PYTHONPATH='/usr/lib/python3/dist-packages',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1'),
                stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        save(B/'launch_receipt.json',dict(controller_pid=p.pid,launched_utc=base.now(),cases=6))
        print(json.dumps(load(B/'launch_receipt.json')))
    elif sys.argv[1:]==['--run']:
        try:run()
        except BaseException:state('CONTROLLER_FAILED',traceback=traceback.format_exc());raise
    else:raise ValueError('--prepare APPROVAL, --launch, or --run')
