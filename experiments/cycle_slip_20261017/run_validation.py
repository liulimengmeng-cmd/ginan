"""Readiness audit. Always preserves failing tests and never starts an AR job."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time


def digest(path):
    value=hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024),b''):
            value.update(chunk)
    return value.hexdigest()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--build',type=Path,required=True)
    parser.add_argument('--evidence',type=Path,required=True)
    parser.add_argument('--checkpoint',type=Path,required=True)
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[2]
    source_files=[
        'src/cpp/common/zhangCycleSlipEvidence.hpp','src/cpp/common/zhangCycleSlipEvidenceState.hpp',
        'src/cpp/common/zhangCycleSlipDopplerWindow.hpp',
        'src/cpp/common/acsQC.cpp','src/cpp/common/satStat.hpp',
        'src/cpp/pea/zhangReference.cpp','src/cpp/pea/zhangReceiverCheckpoint.cpp',
        'src/cpp/pea/zhangReceiverCheckpoint.hpp','src/cpp/pea/zhangCheckpointRuntime.cpp',
        'src/tests/CMakeLists.txt','src/tests/unit/test_ZhangCycleSlipEvidence.cpp',
        'src/tests/unit/r51_slip_qc_production.cpp','src/tests/run_r51_slip_qc.py',
        'experiments/cycle_slip_20261017/run_validation.py',
        'experiments/cycle_slip_20261017/prepare_replay.py',
        'experiments/cycle_slip_20261017/prepare_without_doppler.py',
        'experiments/cycle_slip_20261017/README.md',
        'experiments/cycle_slip_20261017/validate_replay.py',
        'experiments/cycle_slip_20261017/audit_doppler.py',
        'experiments/cycle_slip_20261017/inventory_doppler.py',
        'experiments/cycle_slip_20261017/freeze_offline.py',
        'experiments/cycle_slip_20261017/launch_offline.sh',
        'experiments/cycle_slip_20261017/start_offline.py',
        'experiments/cycle_slip_20261017/acceptance_gates.json']
    source_before={p:digest(root/p) for p in source_files}
    input_files=[Path('/mnt/d/GINAN_R20/inputData/data/MAL200KEN_R_20241990000_01D_30S_MO.rnx'),
        Path('/mnt/d/GINAN_R20/inputData/products/WUM0MGXFIN_20241990000_01D_05M_ORB.SP3'),
        Path('/mnt/d/GINAN_R20/inputData/products/WUM0MGXFIN_20241990000_01D_30S_CLK.CLK'),args.checkpoint]
    input_hashes={str(p):digest(p) for p in input_files}
    args.evidence.mkdir(parents=True,exist_ok=False)
    # Preserve the actual tested candidate, including failures. Hashes alone
    # cannot reconstruct an uncommitted source tree after the next repair.
    for name in source_files:
        destination=args.evidence/'sources'/name
        destination.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(root/name,destination)
    binary=root/'bin/pea';before=digest(binary)
    shutil.copy2(binary,args.evidence/'pea')
    tested_binaries=[root/'bin/zhang_cycle_slip_evidence_tests',root/'bin/zhang_full_rank_tests',
        args.build/'slip_qc_contract/slip_qc_contract']
    tested_binary_hashes={str(p):digest(p) for p in tested_binaries}
    (args.evidence/'test_binaries').mkdir()
    for path in tested_binaries:shutil.copy2(path,args.evidence/'test_binaries'/path.name)
    environment={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1',
                 'ZHANG_R51_RATIO_ONLY':'1','ZHANG_R51_RATIO_THRESHOLD':'2'}
    for flag in ['ZHANG_SLIP_TEST_DOPPLER','ZHANG_SLIP_TEST_IONO_STRESS',
                 'ZHANG_SLIP_TEST_DOPPLER_BLUNDER','ZHANG_SLIP_TEST_PHASE_TREND','ZHANG_SLIP_TEST_TRACE',
                 'ZHANG_SLIP_TEST_NO_DOPPLER']:
        environment.pop(flag,None)
    checks=[]

    def run(name,command,filtered=False,extra_environment=None):
        started=time.monotonic();log=args.evidence/(name+'.log')
        with log.open('w',buffering=1) as output:
            process=subprocess.Popen([str(p) for p in command],stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,text=True,errors='replace',env={**environment,**(extra_environment or {})})
            for line in process.stdout:
                if not filtered or any(token in line for token in
                    ['ZHANG_SLIP','raw_production','PASS','FAIL','failed','missing','incomplete',
                     'false alarm','mismatch','DGEMV','illegal','FATAL','NaN','terminate']):
                    output.write(line)
            code=process.wait()
        checks.append({'name':name,'exit_code':code,'seconds':time.monotonic()-started,'log':str(log)})
        print(json.dumps(checks[-1]),flush=True)
        return code

    # Each Monte Carlo process is single-threaded and owns its detector states.
    # Independent controls can run concurrently without sharing observations.
    cases=[
        ('injected_independent_holdout','--monte-carlo-slip',2000,300000,{'ZHANG_SLIP_TEST_DOPPLER':'1'}),
        ('null_independent_holdout','--monte-carlo',2000,300000,{'ZHANG_SLIP_TEST_DOPPLER':'1'}),
        ('storm_injected_holdout','--monte-carlo-slip',1000,310000,
            {'ZHANG_SLIP_TEST_DOPPLER':'1','ZHANG_SLIP_TEST_IONO_STRESS':'1'}),
        ('storm_null_paired_control','--monte-carlo',1000,310000,
            {'ZHANG_SLIP_TEST_DOPPLER':'1','ZHANG_SLIP_TEST_IONO_STRESS':'1'}),
        ('without_doppler_holdout','--monte-carlo-slip',2000,320000,{}),
        ('without_doppler_paired_control','--monte-carlo',2000,320000,{}),
        ('smooth_phase_trend_holdout','--monte-carlo',1000,330000,
            {'ZHANG_SLIP_TEST_DOPPLER':'1','ZHANG_SLIP_TEST_PHASE_TREND':'1'}),
        ('no_doppler_storm_null','--monte-carlo',1000,340000,{'ZHANG_SLIP_TEST_IONO_STRESS':'1'}),
        ('no_doppler_phase_trend_null','--monte-carlo',1000,350000,{'ZHANG_SLIP_TEST_PHASE_TREND':'1'})]
    with ThreadPoolExecutor(max_workers=min(8,max(2,os.cpu_count() or 2))) as pool:
        jobs=[pool.submit(run,'ctest',['ctest','--test-dir',args.build,
            '--output-on-failure','-j2','--timeout','900'])]
        for name,mode,seeds,offset,flags in cases:
            jobs.append(pool.submit(run,name,[root/'bin/zhang_cycle_slip_evidence_tests',
                mode,str(seeds),str(offset)],extra_environment=flags))
        for job in jobs:job.result()
    fixture=args.evidence/'raw_derived_grid'
    run('network_doppler_inventory',[sys.executable,root/'experiments/cycle_slip_20261017/inventory_doppler.py',
        '--config','/home/rx/GINAN/frozen-r51-satellite-fix-20260918/config/zhang_global_2024199_180_inputs.yaml',
        '--input-root','/mnt/d/GINAN_R20/inputData/data','--output',args.evidence/'network_doppler_inventory.json'])
    run('prepare_raw_fixture',[sys.executable,root/'experiments/cycle_slip_20261017/prepare_replay.py',
        '--input-root','/mnt/d/GINAN_R20/inputData','--output',fixture])
    for dual,nod in [(False,False),(True,False),(False,True),(True,True)]:
        run('raw_derived_grid'+('_dual' if dual else '')+('_no_doppler' if nod else ''),[sys.executable,
            root/'experiments/cycle_slip_20261017/validate_replay.py',
            '--binary',root/'bin/zhang_cycle_slip_evidence_tests','--root',fixture,
            '--grid',*(['--dual'] if dual else []),*(['--without-doppler'] if nod else [])])
    run('network_nis_skip',[root/'bin/zhang_full_rank_tests',
        '--run_test=r51_ratio_policy_retains_structural_rejections','--log_level=test_suite'])
    driver=args.build/'slip_qc_contract/slip_qc_contract'
    run('old_checkpoint_decode',[driver,'--old-checkpoint',args.checkpoint])
    for case in [(0,0,0),(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(1,1,0),(-1,-1,0),
                 (1,1,1),(-1,-1,-1),(9,7,0),(77,60,0),(60,77,0),(0,0,1)]:
        name='raw_MAL2_'+('_'.join(str(n) for n in case))
        run(name,[driver,'--raw','/mnt/d/GINAN_R20/inputData',*case],filtered=True)
    nod_root=args.evidence/'raw_without_doppler'
    run('prepare_production_no_doppler',[sys.executable,root/'experiments/cycle_slip_20261017/prepare_without_doppler.py',
        '--input-root','/mnt/d/GINAN_R20/inputData','--output-root',nod_root])
    for case in [(0,0,0),(1,0,0),(0,1,0),(1,1,0),(1,1,1),(-1,-1,-1),(0,0,1)]:
        name='raw_noD_MAL2_'+('_'.join(str(n) for n in case))
        run(name,[driver,'--raw',nod_root,*case],filtered=True)
    regressions=[
        (110314,60,{'ZHANG_SLIP_TEST_DOPPLER':'1'}),(110494,60,{'ZHANG_SLIP_TEST_DOPPLER':'1'}),
        (160402,60,{'ZHANG_SLIP_TEST_DOPPLER':'1'}),(160482,60,{'ZHANG_SLIP_TEST_DOPPLER':'1'}),
        (130866,120,{'ZHANG_SLIP_TEST_DOPPLER':'1','ZHANG_SLIP_TEST_IONO_STRESS':'1'}),
        (170418,120,{'ZHANG_SLIP_TEST_DOPPLER':'1','ZHANG_SLIP_TEST_IONO_STRESS':'1'}),
        (140978,330,{}),(180234,330,{}),(180270,330,{}),(180762,330,{}),(180886,330,{})]
    for seed,limit,flags in regressions:
        run('long_tail_'+str(seed),[root/'bin/zhang_cycle_slip_evidence_tests',
            '--monte-carlo-slip','1',str(seed)],extra_environment=flags)
    source_unchanged=source_before=={p:digest(root/p) for p in source_files}
    binary_unchanged=digest(binary)==before
    tested_binaries_unchanged=tested_binary_hashes=={str(p):digest(p) for p in tested_binaries}
    snapshot_matches=all(digest(args.evidence/'sources'/p)==h for p,h in source_before.items())
    inputs_unchanged=input_hashes=={str(p):digest(p) for p in input_files}
    nod_logs=list(args.evidence.glob('raw_noD_MAL2_*.log'))
    nod_trace_consistent=len(nod_logs)==7 and all('ZHANG_SLIP_EVIDENCE' in p.read_text() and
        not re.search(r'ZHANG_SLIP_EVIDENCE[^\n]*\bdoppler=1\b',p.read_text()) for p in nod_logs)
    def summary(name):
        try:return json.loads((args.evidence/(name+'.log')).read_text().splitlines()[-1])
        except (ValueError,IndexError,FileNotFoundError):return None
    performance=[]
    for name,minimum_current,maximum_delay in [('injected_independent_holdout',.99,60),
            ('storm_injected_holdout',.98,120),('without_doppler_holdout',.85,330)]:
        data=summary(name)
        control_name={'injected_independent_holdout':'null_independent_holdout',
            'storm_injected_holdout':'storm_null_paired_control',
            'without_doppler_holdout':'without_doppler_paired_control'}[name]
        control=summary(control_name)
        if not data or not control:
            performance.append({'name':name,'passed':False,'reason':'MISSING_MACHINE_READABLE_SUMMARY'});continue
        already=set(control['event_epoch_excluded_seed_offsets'])
        newly=set(data['event_epoch_excluded_seed_offsets'])-already
        rate=len(newly)/max(1,data['injected_cases']-len(already))
        performance.append({'name':name,'causal_current_epoch_rate':rate,
            'matched_null_excluded_at_event':len(already),'new_current_epoch_exclusions':len(newly),
            'maximum_delay_seconds':data['max_detection_delay_s'],'required_current_rate':minimum_current,
            'allowed_maximum_delay_seconds':maximum_delay,
            'additional_retirement_events':data['retirements']-data['retired_cases'],
            'passed':rate>=minimum_current and data['max_detection_delay_s']<=maximum_delay and
                data['retirements']==data['retired_cases']})
    long_tail=[]
    for seed,limit,flags in regressions:
        data=summary('long_tail_'+str(seed))
        long_tail.append({'seed':seed,'allowed_maximum_delay_seconds':limit,
            'maximum_delay_seconds':data['max_detection_delay_s'] if data else None,
            'passed':bool(data and data['detected_cases']==1 and data['already_excluded_before_injection']==0 and
                0<=data['max_detection_delay_s']<=limit and data['wrong_integer_certificates']==0 and data['wrong_retirement_masks']==0)})
    result={'checks':checks,'binary_sha256':before,'binary_unchanged_during_tests':digest(binary)==before,
        'candidate_snapshot':str(args.evidence/'sources'),
        'candidate_snapshot_hashes_match':snapshot_matches,
        'tested_binary_hashes':tested_binary_hashes,'tested_binaries_unchanged':tested_binaries_unchanged,
        'source_hashes':source_before,'source_unchanged_during_tests':source_unchanged,
        'checkpoint_sha256':digest(args.checkpoint),
        'input_hashes':input_hashes,'inputs_unchanged':inputs_unchanged,
        'no_doppler_production_trace_consistent':nod_trace_consistent,
        'long_tail_regression_gates':long_tail,
        'performance_gates':performance,
        'all_execution_checks_passed':all(c['exit_code']==0 for c in checks) and source_unchanged and
            binary_unchanged and tested_binaries_unchanged and snapshot_matches and inputs_unchanged and nod_trace_consistent and
            all(p['passed'] for p in performance) and all(p['passed'] for p in long_tail),
        'ar_started':False,'restart_approved':False,
        'limitations':['No full network KF/graph/product replay claimed by the QC harness',
            'Raw LLI=0 does not establish perfect no-slip truth',
            'No-slip Gaussian and injected-noise rates must be reported separately',
            'Current-epoch detection is not guaranteed for every noise realization or absent-Doppler link',
            '121 of 180 input headers declare both selected GPS Doppler fields; per-epoch availability can be lower',
            'Quarantine is not successful identification or correction of a cycle slip']}
    (args.evidence/'validation_result.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({'all_execution_checks_passed':result['all_execution_checks_passed'],
        'ar_started':False,'report':str(args.evidence/'validation_result.json')}),flush=True)
    return 0 if result['all_execution_checks_passed'] else 1


if __name__=='__main__':raise SystemExit(main())
