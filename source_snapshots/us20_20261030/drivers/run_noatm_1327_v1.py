"""Isolated Galileo replay: disable only regional ionosphere and ZTD constraints."""
from pathlib import Path
import copy, importlib.util, io, json, os, shutil, subprocess, sys, tarfile, traceback

OLD=Path('/home/rx/GINAN/us20_osb_event_1327_galileo_ratio3_1h_20261029_v1')
BASE=Path('/home/rx/GINAN/us20_osb_event_1327_1h_20261029_v1')
B=Path('/home/rx/GINAN/us20_osb_event_1327_galileo_noatm_ratio3_1h_20261030_v1')
spec=importlib.util.spec_from_file_location('frozen_one_hour_harness',BASE/'scripts/run_osb1327_one_hour_v1.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);m.B=B
load,save,sha=m.load,m.save,m.sha

OLD_EXPECTATION="""    assert applied
    applied_by_system=collections.Counter(e.get('system','UNSPECIFIED') for e in applied)
    assert set(applied_by_system)=={'GAL'},applied_by_system
    assert all(e['full_joint_R_shared']==e['receiver_common_GF_direction_preserved']=='1' and e['user_truth_geometry_used']=='0' for e in applied)
"""
NEW_EXPECTATION="""    assert not applied, 'ATMOSPHERIC_CONSTRAINTS_MUST_BE_DISABLED'
    applied_by_system=collections.Counter(e.get('system','UNSPECIFIED') for e in applied)
    assert not applied_by_system
"""
OLD_CAPTURE="    assert all(e['marker_matched']==e['full_CP_PHASE_IONO_R']=='1' and e['actual_rows']==e['expected_rows'] for e in captures)"
NEW_CAPTURE=OLD_CAPTURE+"\n    assert all(e['actual_rows']==e['expected_rows']==e['actual_trop_rows']==e['expected_trop_rows']=='0' for e in captures), 'ATMOSPHERIC_CAPTURE_ROWS_MUST_BE_ZERO'"

def bind_validator():
    spec=importlib.util.spec_from_file_location('noatm_calendar_validator',B/'scripts/calendar_accounting_galileo_noatm_v1.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);m.accounting=module

def prepare():
    assert not B.exists(),'Inspect existing state/PIDs before resuming; never duplicate runs'
    prior=load(OLD/'experiment.json')
    assert prior['ratio_threshold']==3 and prior['observation_systems']==['GAL']
    assert load(OLD/'analysis/summary.json')['completed']
    assert sha(prior['binary'])==prior['binary_sha256']=='fe75641d62a7abb9477f0c7d2b541e491ebc209b5f6dad7817ebdb454d90f9a8'
    assert len(os.sched_getaffinity(0))>=6 and shutil.disk_usage(B.parent).free>25*1024**3
    B.mkdir();(B/'scripts').mkdir()
    shutil.copy2(__file__,B/'scripts/run_noatm_1327_v1.py')
    shutil.copy2(BASE/'scripts/run_osb1327_one_hour_v1.py',B/'scripts/run_osb1327_one_hour_v1.py')
    shutil.copy2(OLD/'reference_evaluation_only.json',B/'reference_evaluation_only.json')
    v0=OLD/'scripts/calendar_accounting_galileo_ratio3_v1.py';source=v0.read_text()
    assert source.count(OLD_EXPECTATION)==source.count(OLD_CAPTURE)==1
    derivative=B/'scripts/calendar_accounting_galileo_noatm_v1.py'
    derivative.write_text(source.replace(OLD_EXPECTATION,NEW_EXPECTATION).replace(OLD_CAPTURE,NEW_CAPTURE))
    restored=derivative.read_text().replace(NEW_CAPTURE,OLD_CAPTURE).replace(NEW_EXPECTATION,OLD_EXPECTATION)
    assert restored==source
    save(B/'validator_derivation.json',dict(original_sha256=sha(v0),isolated_sha256=sha(derivative),
        changed_expectations=['Regional ionosphere applied count must be zero','Ionosphere/ZTD captured pseudo rows must be zero'],
        all_other_calendar_commit_FLOAT_authority_ratio_POS_requirements_unchanged=True))
    save(B/'user_approval.json',dict(request='重新进行单系统伽利略试验；关闭对流层和电离层约束，绘制两策略各自AR/FLOAT及两策略FLOAT三组ENU图',
        scope='Same three stations, two history policies, 13:10--14:10 GPST, ratio>3, 7.5 deg',
        regional_IONO='OFF',regional_TROP='OFF',physical_models_and_STEC_ZTD_estimation_retained=True,
        server_products_unchanged=True,provider_rerun=False,created_utc=m.base.now()))
    rows=[];diffs=[]
    for old in prior['cases']:
        for name,digest in [('config','config_sha256'),('environment','environment_sha256'),('observation','observation_sha256')]:
            assert sha(old[name])==old[digest],(old['station'],name)
        c0=m.yaml.safe_load(Path(old['config']).read_text());e0=load(old['environment'])
        c=copy.deepcopy(c0);e=copy.deepcopy(e0)
        assert e0['ZHANG_USER_REGIONAL_IONO']==e0['ZHANG_USER_REGIONAL_TROP']=='ON'
        assert e0['ZHANG_USER_REGIONAL_MODE']=='SD_PLANE_JOINT'
        e['ZHANG_USER_REGIONAL_IONO']=e['ZHANG_USER_REGIONAL_TROP']='OFF'
        out=B/'cases'/old['policy']/old['case'];out.mkdir(parents=True)
        c['outputs']['outputs_root']=str(out)
        c['outputs']['metadata']['config_description']+='_NO_ATM_CONSTRAINTS'
        z=c['processing_options']['gnss_general']['zhang_pppar']
        assert z['checkpoint_restore_path']=='' and not z['deterministic_checkpoint']
        z['checkpoint_output_directory']=str(out/'checkpoints');z['checkpoint_runtime_id']+='_NO_ATM'
        e['ZHANG_RESET_AUDIT_DIRECTORY']=str(out/'state_audit')
        e['ZHANG_USER_TROP_CONSUMPTION_LEDGER']=str(out/'trop_consumption_ledger')
        assert c['processing_options']['ambiguity_resolution']['solution_ratio_threshold']==3
        assert c['processing_options']['ambiguity_resolution']['elevation_mask']==7.5
        assert {k for k,v in c['processing_options']['gnss_general']['sys_options'].items() if v['process']}=={'gal'}
        assert c['estimation_parameters']['receivers']['global']['ion_stec']['estimated']==[True]
        assert c['estimation_parameters']['receivers']['global']['trop']['estimated']==[True]
        assert c['processing_options']['ppp_filter']['ionospheric_components']['corr_mode']=='estimate'
        assert c['receiver_options']['global']['models']['troposphere']['enable']
        assert c['receiver_options']['global']['models']['ionospheric_components']['enable']
        restored=copy.deepcopy(c);restored['outputs']=copy.deepcopy(c0['outputs'])
        for key in ['checkpoint_output_directory','checkpoint_runtime_id']:
            restored['processing_options']['gnss_general']['zhang_pppar'][key]=c0['processing_options']['gnss_general']['zhang_pppar'][key]
        assert restored==c0
        restored_env=copy.deepcopy(e)
        changed_keys=['ZHANG_RESET_AUDIT_DIRECTORY','ZHANG_USER_TROP_CONSUMPTION_LEDGER','ZHANG_USER_REGIONAL_IONO','ZHANG_USER_REGIONAL_TROP']
        for k in changed_keys:restored_env[k]=e0[k]
        assert restored_env==e0
        (out/'execution.yaml').write_text(m.yaml.safe_dump(c,sort_keys=False));save(out/'environment.json',e)
        item=copy.deepcopy(old)
        for key in list(item):
            if key.startswith('scientific_settings_equal_except_'):item.pop(key)
        item.update(working_directory=str(out),config=str(out/'execution.yaml'),config_sha256=sha(out/'execution.yaml'),
            environment=str(out/'environment.json'),environment_sha256=sha(out/'environment.json'),
            command=[prior['binary'],'--rnx_inputs',old['observation'],'-y',str(out/'execution.yaml')],
            original_config=old['config'],original_config_sha256=old['config_sha256'],
            original_environment=old['environment'],original_environment_sha256=old['environment_sha256'],
            scientific_settings_equal_except_regional_IONO_TROP_switches_and_output=True)
        rows.append(item)
        save(out/'runtime.json',dict(state='QUEUED',case=item['case'],station=item['station'],policy=item['policy']))
        diffs.append(dict(station=item['station'],policy=item['policy'],only_scientific_environment_changes={
            'ZHANG_USER_REGIONAL_IONO':['ON','OFF'],'ZHANG_USER_REGIONAL_TROP':['ON','OFF']},
            execution_model_unchanged=True,unchanged_joint_clock_phase_covariance_path=True))
    exp=copy.deepcopy(prior)
    exp.update(schema='US20_OSB1327_GALILEO_NO_ATMOSPHERIC_CONSTRAINTS_RATIO3_V1',created_utc=m.base.now(),cases=rows,
        baseline_root=str(OLD),baseline_experiment_sha256=sha(OLD/'experiment.json'),
        approval_sha256=sha(B/'user_approval.json'),regional_ionosphere_constraint=False,
        regional_troposphere_constraint=False,physical_STEC_ZTD_estimation=True,
        provider_rerun=False,frozen_harness_sha256=sha(BASE/'scripts/run_osb1327_one_hour_v1.py'))
    save(B/'experiment.json',exp)
    save(B/'noatm_derivation.json',dict(passed=True,cases=diffs,binary_sha256=exp['binary_sha256'],
        reference_sha256=sha(B/'reference_evaluation_only.json'),source_code_modified=False,
        physical_atmosphere_models_and_stochastic_initialisation_unchanged=True))
    bind_validator()
    # A genuine atmosphere-constrained old case must fail the new zero-row expectation.
    old=prior['cases'][0]
    try:m.accounting.validate_case_calendar(Path(old['working_directory']),old,m.START,120)
    except AssertionError:
        failure=traceback.format_exc();assert 'ATMOSPHERIC_CONSTRAINTS_MUST_BE_DISABLED' in failure
    else:raise AssertionError('Constrained baseline incorrectly accepted as no-atmospheric-constraint experiment')
    save(B/'validator_negative_control.json',dict(passed=True,constrained_baseline_rejected=True,reason=failure))
    m.state('PREPARED_NO_ATMOSPHERIC_CONSTRAINTS',cases=6,new_solver_processes=0)
    print(json.dumps(dict(prepared=True,root=str(B),cases=6,IONO='OFF',TROP='OFF',binary_unchanged=True)),flush=True)

def native_audit(exp):
    audit=[]
    for item in exp['cases']:
        out=Path(item['working_directory']);ev=load(out/'event_receipts.json')
        ppp=ev['ZHANG_USER_PPP_EPOCH'];pol=ev['ZHANG_USER_GONSIS_POLICY'];ar=ev['ZHANG_USER_GONSIS_AR_SUMMARY']
        assert len(ppp)==120 and len({x['time'] for x in ppp})==120
        assert all(x['ratio_threshold']=='3' and x['ratio_comparison']=='STRICT_GREATER' for x in pol)
        assert {x['time'] for x in pol}=={x['time'] for x in ppp if x['ar_authorized']=='1'}
        assert not ev['ZHANG_USER_REGIONAL_EPOCH']
        caps=ev['ZHANG_USER_GONSIS_REGIONAL_CAPTURE']
        assert all(x['actual_rows']==x['expected_rows']==x['actual_trop_rows']==x['expected_trop_rows']=='0' for x in caps)
        assert all(x['system']=='GAL' or int(x['eligible_dimension'])==int(x['final_rank'])==0 for x in ar)
        assert all(x['satellite'].startswith('E') for x in ev['ZHANG_USER_SHARED_ENDPOINT'] if x.get('search_eligible')=='1')
        forbidden=['ZHANG_USER_REGIONAL_SD ','ZHANG_USER_REGIONAL_TROP ','ZHANG_USER_REGIONAL_EPOCH ']
        raw_counts={tag:0 for tag in forbidden};postfix_counts=[]
        for p in out.glob('Network-*.TRACE'):
            for line in p.open():
                for tag in forbidden:
                    if line.startswith(tag):raw_counts[tag]+=1
                if line.startswith('ZHANG_USER_GONSIS_OUTPUT_RESIDUAL '):
                    e=m.parser.event(line)
                    assert e['regional_iono_rows']==e['regional_trop_rows']=='0'
                    postfix_counts.append(dict(time=e['time'],iono_rows=0,trop_rows=0))
        assert not any(raw_counts.values())
        ledger=out/'trop_consumption_ledger'
        assert not ledger.exists() or not any(ledger.rglob('*'))
        snapshots=sorted((out/'state_audit').glob('*/committed_posterior/state_keys.txt'))
        fresh={x['time'] for x in ppp if x['status']=='MEASUREMENT_UPDATED'}
        snapshot_times=set();state_counts=[]
        for p in snapshots:
            keys=m.keys(p);assert all(not k[2] or k[2].startswith('E') for k in keys)
            ztd=[k for k in keys if k[0]==21]
            ion=[k for k in keys if k[0]==25]
            assert len(ztd)>=1 and len(ion)>=1
            snapshot_times.add(str(m.dt.datetime(1980,1,6)+m.dt.timedelta(seconds=int(p.parents[1].name))))
            state_counts.append(dict(epoch_seconds=p.parents[1].name,trop_states=len(ztd),STEC_states=len(ion),state_key_types=sorted(set(k[0] for k in keys))))
        assert fresh<=snapshot_times
        audit.append(dict(station=item['station'],policy=item['policy'],native_PPP_epochs=120,
            fresh_updates=len(fresh),authorized_AR_epochs=len(pol),posterior_snapshots=len(snapshots),
            non_Galileo_satellite_keys=0,regional_iono_applied_rows=0,regional_trop_applied_rows=0,
            atmospheric_trace_counts=raw_counts,postfix_output_checks=len(postfix_counts),
            zero_pseudo_row_captures=len(caps),state_counts=state_counts))
    return dict(passed=True,cases=audit,regional_IONO='OFF',regional_TROP='OFF',service_products_remain_GEC=True)

def pack(exp,audit):
    save(B/'noatm_native_system_audit.json',audit)
    # Keep the old plot loader's clear Galileo-only posterior assertion interface.
    save(B/'galileo_native_system_audit.json',dict(passed=True,cases=audit['cases'],no_atmospheric_constraints=True))
    files={}
    with tarfile.open(B/'delivery_results.tar.gz','r:gz') as tf:
        old=json.load(tf.extractfile('payload_manifest.json'))['files']
        files.update({n:B/n for n in old})
    for n in ['noatm_derivation.json','noatm_native_system_audit.json','galileo_native_system_audit.json',
              'validator_derivation.json','validator_negative_control.json','scripts/run_noatm_1327_v1.py',
              'scripts/calendar_accounting_galileo_noatm_v1.py']:
        files[n]=B/n
    archive=B/'noatm_results.tar.gz';assert not archive.exists();payload={}
    with tarfile.open(archive,'w:gz',compresslevel=3) as tf:
        for n,p in files.items():
            data=p.read_bytes();payload[n]=dict(bytes=len(data),sha256=sha(p));i=tarfile.TarInfo(n);i.size=len(data);tf.addfile(i,io.BytesIO(data))
        data=json.dumps(dict(files=payload),indent=2).encode();i=tarfile.TarInfo('payload_manifest.json');i.size=len(data);tf.addfile(i,io.BytesIO(data))
    m.state('COMPLETE_NOATM_REMOTE_AWAITING_LOCAL_DELIVERY',completed_cases=6,archive=str(archive),
        archive_sha256=sha(archive),archive_bytes=archive.stat().st_size,atmospheric_pseudo_rows=0,native_audit_passed=True)

def run():
    # The frozen harness starts only from PREPARED; the launch wrapper performs this transition once.
    bind_validator();m.run()
    assert load(B/'workflow_state.json')['state']=='COMPLETE_REMOTE_ANALYSIS_AWAITING_LOCAL_DELIVERY'
    exp=load(B/'experiment.json');audit=native_audit(exp)
    assert sha(OLD/'experiment.json')==exp['baseline_experiment_sha256']
    for x in exp['cases']:
        assert sha(x['original_config'])==x['original_config_sha256']
        assert sha(x['original_environment'])==x['original_environment_sha256']
    pack(exp,audit)

if __name__=='__main__':
    if sys.argv[1:]==['--prepare']:prepare()
    elif sys.argv[1:]==['--launch']:
        assert load(B/'workflow_state.json')['state']=='PREPARED_NO_ATMOSPHERIC_CONSTRAINTS'
        assert not (B/'launch_receipt.json').exists()
        m.state('PREPARED',cases=6)
        with (B/'controller.log').open('x') as log:
            p=subprocess.Popen([m.PY,str(B/'scripts/run_noatm_1327_v1.py'),'--run'],cwd=B,
                env=dict(os.environ,PYTHONPATH='/usr/lib/python3/dist-packages',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1'),
                stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        save(B/'launch_receipt.json',dict(controller_pid=p.pid,launched_utc=m.base.now(),cases=6,IONO='OFF',TROP='OFF'))
        print(json.dumps(load(B/'launch_receipt.json')),flush=True)
    elif sys.argv[1:]==['--run']:
        try:run()
        except BaseException:m.state('NOATM_CONTROLLER_FAILED',traceback=traceback.format_exc());raise
    else:raise ValueError('--prepare, --launch, or --run')
