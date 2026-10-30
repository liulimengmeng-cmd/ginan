"""Galileo-only user replay against unchanged 20-station GEC products."""
from pathlib import Path
import copy,datetime as dt,importlib.util,io,json,os,shutil,subprocess,sys,tarfile,traceback
OLD=Path('/home/rx/GINAN/us20_osb_event_1327_ratio3_1h_20261029_v1')
BASE=Path('/home/rx/GINAN/us20_osb_event_1327_1h_20261029_v1')
B=Path('/home/rx/GINAN/us20_osb_event_1327_galileo_ratio3_1h_20261029_v1')
spec=importlib.util.spec_from_file_location('frozen_window_harness',BASE/'scripts/run_osb1327_one_hour_v1.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);m.B=B
load,save,sha=m.load,m.save,m.sha

def bind_validator():
    spec=importlib.util.spec_from_file_location('galileo_calendar_validator',B/'scripts/calendar_accounting_galileo_ratio3_v1.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);m.accounting=module

def prepare():
    assert not B.exists(),'Inspect existing workflow and actual PIDs before resuming'
    previous=load(OLD/'experiment.json');assert previous['ratio_threshold']==3
    assert load(OLD/'analysis/summary.json')['completed']
    assert len(os.sched_getaffinity(0))>=6 and shutil.disk_usage(B.parent).free>25*1024**3
    assert sha(previous['binary'])==previous['binary_sha256']
    B.mkdir();(B/'scripts').mkdir()
    shutil.copy2(__file__,B/'scripts/run_galileo_only_1327_v1.py')
    shutil.copy2(BASE/'scripts/run_osb1327_one_hour_v1.py',B/'scripts/run_osb1327_one_hour_v1.py')
    shutil.copy2(OLD/'reference_evaluation_only.json',B/'reference_evaluation_only.json')
    validator=OLD/'scripts/calendar_accounting_ratio3_v1.py';source=validator.read_text()
    original="assert set(applied_by_system)=={'GPS','GAL','BDS'},applied_by_system"
    replacement="assert set(applied_by_system)=={'GAL'},applied_by_system"
    assert source.count(original)==1
    derivative=B/'scripts/calendar_accounting_galileo_ratio3_v1.py'
    derivative.write_text(source.replace(original,replacement))
    assert derivative.read_text().replace(replacement,original)==source
    assert "e['ratio_threshold']=='3'" in derivative.read_text()
    save(B/'validator_derivation.json',dict(original_sha256=sha(validator),isolated_sha256=sha(derivative),
        only_changed_expectation='Applied regional systems must equal GAL instead of GPS/GAL/BDS',
        ratio_requirement_3_unchanged=True,all_remaining_acceptance_source_unchanged=True))
    save(B/'user_approval.json',dict(request='请你绘制enu图后，再进行实验，只用伽利略去进行用户端的测试',
        scope='Same three stations and two history policies, 13:10--14:10 GPST; ratio>3; Galileo-only USER observations',
        provider_remains_GEC=True,ENU_figures_completed_before_solver_launch=True,created_utc=m.base.now()))
    rows=[];diffs=[]
    for old in previous['cases']:
        assert sha(old['config'])==old['config_sha256'] and sha(old['environment'])==old['environment_sha256']
        assert sha(old['observation'])==old['observation_sha256']
        cfg0=m.yaml.safe_load(Path(old['config']).read_text());env0=load(old['environment'])
        cfg=copy.deepcopy(cfg0);env=copy.deepcopy(env0)
        out=B/'cases'/old['policy']/old['case'];out.mkdir(parents=True)
        options=cfg['processing_options']['gnss_general']['sys_options']
        assert {k for k,v in options.items() if v['process']}=={'gps','gal','bds'}
        assert cfg['processing_options']['ambiguity_resolution']['solution_ratio_threshold']==3
        for system,v in options.items():
            v['process']=v['ambiguity_resolution']=(system=='gal')
            if system!='gal':
                for key in ['receiver_amb_pivot','network_amb_pivot']:
                    if key in v:v[key]=False
        assert options['gal']['code_priorities']==['L1C','L5Q']
        cfg['receiver_options']['global']['rec_reference_system']='GAL'
        assert not any(v.get('rec_reference_system','GAL')!='GAL' for k,v in cfg['receiver_options'].items() if k!='global')
        cfg['outputs']['outputs_root']=str(out)
        cfg['outputs']['metadata']['config_description']+='_GAL_ONLY'
        z=cfg['processing_options']['gnss_general']['zhang_pppar']
        assert z['checkpoint_restore_path']=='' and not z['deterministic_checkpoint']
        z['checkpoint_output_directory']=str(out/'checkpoints');z['checkpoint_runtime_id']+='_GAL_ONLY'
        env['ZHANG_RESET_AUDIT_DIRECTORY']=str(out/'state_audit')
        env['ZHANG_USER_TROP_CONSUMPTION_LEDGER']=str(out/'trop_consumption_ledger')
        restored=copy.deepcopy(cfg);restored['outputs']=copy.deepcopy(cfg0['outputs'])
        restored['processing_options']['gnss_general']['sys_options']=copy.deepcopy(cfg0['processing_options']['gnss_general']['sys_options'])
        restored['receiver_options']['global']['rec_reference_system']=cfg0['receiver_options']['global']['rec_reference_system']
        for key in ['checkpoint_output_directory','checkpoint_runtime_id']:
            restored['processing_options']['gnss_general']['zhang_pppar'][key]=cfg0['processing_options']['gnss_general']['zhang_pppar'][key]
        assert restored==cfg0
        restored_env=copy.deepcopy(env)
        for key in ['ZHANG_RESET_AUDIT_DIRECTORY','ZHANG_USER_TROP_CONSUMPTION_LEDGER']:restored_env[key]=env0[key]
        assert restored_env==env0
        (out/'execution.yaml').write_text(m.yaml.safe_dump(cfg,sort_keys=False));save(out/'environment.json',env)
        item=copy.deepcopy(old)
        item.update(working_directory=str(out),config=str(out/'execution.yaml'),config_sha256=sha(out/'execution.yaml'),
            environment=str(out/'environment.json'),environment_sha256=sha(out/'environment.json'),
            command=[previous['binary'],'--rnx_inputs',old['observation'],'-y',str(out/'execution.yaml')],
            original_config=old['config'],original_config_sha256=old['config_sha256'],original_environment_sha256=old['environment_sha256'],
            scientific_settings_equal_except_constellation_clock_reference_and_output=True)
        item.pop('scientific_settings_equal_except_ratio_and_output',None);rows.append(item)
        save(out/'runtime.json',dict(state='QUEUED',case=item['case'],station=item['station'],policy=item['policy']))
        diffs.append(dict(station=item['station'],policy=item['policy'],original_config_sha256=old['config_sha256'],
            new_config_sha256=item['config_sha256'],enabled_observation_systems=['GAL'],enabled_AR_systems=['GAL'],
            receiver_clock_reference_system='GAL',only_system_mask_clock_reference_and_output_paths_changed=True))
    exp=copy.deepcopy(previous)
    exp.update(schema='US20_OSB1327_GALILEO_ONLY_RATIO3_ONE_HOUR_V1',created_utc=m.base.now(),cases=rows,
        baseline_root=str(OLD),baseline_experiment_sha256=sha(OLD/'experiment.json'),approval_sha256=sha(B/'user_approval.json'),
        observation_systems=['GAL'],user_AR_systems=['GAL'],receiver_reference_system='GAL',
        provider_product_systems=['GPS','GAL','BDS'],provider_rerun=False,ratio_threshold=3)
    save(B/'experiment.json',exp)
    save(B/'galileo_derivation.json',dict(passed=True,case_differences=diffs,unchanged_binary_sha256=exp['binary_sha256'],
        unchanged_reference_sha256=sha(B/'reference_evaluation_only.json'),unchanged_provider_export_sha256=exp['provider_export_sha256'],
        source_model_modified=False,ratio_and_all_other_numeric_thresholds_unchanged=True))
    bind_validator()
    # The genuine GEC control must be rejected by a Galileo-only regional-system audit.
    old=previous['cases'][0]
    try:m.accounting.validate_case_calendar(Path(old['working_directory']),old,m.START,120)
    except AssertionError:
        failure=traceback.format_exc();assert replacement in failure
    else:raise AssertionError('GEC baseline incorrectly passed Galileo-only system audit')
    save(B/'validator_negative_control.json',dict(passed=True,GEC_control_rejected=True,reason=failure,new_solver_processes=0))
    m.state('PREPARED',cases=6,observation_systems=['GAL'],ratio_threshold=3)
    print(json.dumps(dict(prepared=True,root=str(B),cases=6,observation_systems=['GAL'],receiver_clock_reference='GAL',ratio=3)),flush=True)

def run():
    bind_validator();m.run()
    assert load(B/'workflow_state.json')['state']=='COMPLETE_REMOTE_ANALYSIS_AWAITING_LOCAL_DELIVERY'
    exp=load(B/'experiment.json');native=[]
    for item in exp['cases']:
        root=Path(item['working_directory']);ev=load(root/'event_receipts.json');ppp=ev['ZHANG_USER_PPP_EPOCH'];pol=ev['ZHANG_USER_GONSIS_POLICY']
        assert len(ppp)==120
        assert len(pol)==len({r['time'] for r in pol})
        assert {r['time'] for r in pol}=={r['time'] for r in ppp if r['ar_authorized']=='1'}
        assert all(r['ratio_threshold']=='3' and r['ratio_comparison']=='STRICT_GREATER' for r in pol)
        regional={r['system'] for r in ev['ZHANG_USER_REGIONAL_EPOCH'] if r.get('status')=='APPLIED'}
        assert regional=={'GAL'}
        assert all(r['system']=='GAL' or int(r['eligible_dimension'])==int(r['final_rank'])==0 for r in ev['ZHANG_USER_GONSIS_AR_SUMMARY'])
        assert all(r['satellite'].startswith('E') for r in ev['ZHANG_USER_SHARED_ENDPOINT'] if r.get('search_eligible')=='1')
        snapshots=sorted((root/'state_audit').glob('*/committed_posterior/state_keys.txt'));assert len(snapshots)==120
        for path in snapshots:
            keys=m.keys(path);assert all(not k[2] or k[2].startswith('E') for k in keys),path
        native.append(dict(station=item['station'],policy=item['policy'],native_PPP_epochs=120,
            fresh_updates=load(root/'runtime_acceptance.json')['fresh_updates'],authorized_AR_epochs=len(pol),
            native_ratio_threshold=3,regional_applied_systems=['GAL'],posterior_snapshots_Galileo_only=120,
            non_Galileo_satellite_state_keys=0,non_Galileo_eligible_AR_dimensions=0))
        assert sha(item['original_config'])==item['original_config_sha256']
    assert sha(OLD/'experiment.json')==exp['baseline_experiment_sha256']
    save(B/'galileo_native_system_audit.json',dict(passed=True,cases=native,service_products_remain_GEC=True))
    files={}
    with tarfile.open(B/'delivery_results.tar.gz','r:gz') as tar:
        files.update({name:B/name for name in json.load(tar.extractfile('payload_manifest.json'))['files']})
    for name in ['galileo_derivation.json','galileo_native_system_audit.json','validator_derivation.json','validator_negative_control.json',
                 'scripts/run_galileo_only_1327_v1.py','scripts/calendar_accounting_galileo_ratio3_v1.py']:
        files[name]=B/name
    archive=B/'galileo_results.tar.gz';assert not archive.exists();payload={}
    with tarfile.open(archive,'w:gz',compresslevel=3) as tar:
        for name,path in files.items():
            data=path.read_bytes();payload[name]=dict(bytes=len(data),sha256=sha(path));info=tarfile.TarInfo(name);info.size=len(data);tar.addfile(info,io.BytesIO(data))
        data=json.dumps(dict(files=payload),indent=2).encode();info=tarfile.TarInfo('payload_manifest.json');info.size=len(data);tar.addfile(info,io.BytesIO(data))
    m.state('COMPLETE_GALILEO_ONLY_REMOTE_AWAITING_LOCAL_DELIVERY',completed_cases=6,observation_systems=['GAL'],ratio_threshold=3,
        archive=str(archive),archive_sha256=sha(archive),archive_bytes=archive.stat().st_size,native_system_audit_passed=True)

if __name__=='__main__':
    if sys.argv[1:]==['--prepare']:prepare()
    elif sys.argv[1:]==['--launch']:
        assert load(B/'workflow_state.json')['state']=='PREPARED' and not (B/'launch_receipt.json').exists()
        with (B/'controller.log').open('x') as log:
            p=subprocess.Popen([m.PY,str(B/'scripts/run_galileo_only_1327_v1.py'),'--run'],cwd=B,
                env=dict(os.environ,PYTHONPATH='/usr/lib/python3/dist-packages',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1'),
                stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        save(B/'launch_receipt.json',dict(controller_pid=p.pid,launched_utc=m.base.now(),cases=6,observation_systems=['GAL']))
        print(json.dumps(load(B/'launch_receipt.json')),flush=True)
    elif sys.argv[1:]==['--run']:
        try:run()
        except BaseException:m.state('GALILEO_ONLY_CONTROLLER_FAILED',traceback=traceback.format_exc());raise
    else:raise ValueError('--prepare, --launch, or --run')
