"""Audit completed controlled event replay. Never authorizes a run still in progress.

Arguments: evidence-dir build-dir [--completed-hour]. The optional mode audits
the already completed user-authorized hour without authorizing another run.
A failed/missing gate remains false, never NOT_EVALUATED=pass.
"""
from pathlib import Path
import csv
from datetime import datetime,timedelta
import hashlib
import json
import math
import re
import statistics
import subprocess
import sys
import xml.etree.ElementTree as ET

evidence,build=map(lambda s:Path(s).resolve(),sys.argv[1:3])
completed_hour=len(sys.argv)==4 and sys.argv[3]=='--completed-hour'
assert len(sys.argv)==3 or completed_hour
gate_name='completed_hour_gates' if completed_hour else 'formal_gates'
report_name='completed_hour_event_audit' if completed_hour else 'event_audit'
# Re-auditing can never leave an old successful authorization behind on error.
(evidence/(gate_name+'.json')).write_text(json.dumps(dict(all_required_tests_passed=False,
    event_0206_verified=False,event_0237_verified=False,status='AUDIT_PENDING')))
state=json.loads((evidence/'run_state.json').read_text())
run=Path(state['run_root'])
assert (state['status']=='EXITED' and state['returncode']==0) if completed_hour else \
    (state['status']=='CONTROLLED_STOPPED' and state['returncode']==130)
assert not state.get('solver_still_running',True)
report=dict(scope='completed_in_place_hour_event_consistency_audit' if completed_hour else
    'controlled_network_replay_not_completed_formal_hour',
    binary_sha256=state['binary_sha256'],run_root=str(run),errors=[],epochs=[],events={})
def error(message): report['errors'].append(message)

candidate=Path(state['source_binary']).parent
manifest=json.loads((candidate/'manifest.json').read_text())
if hashlib.sha256(Path(state['binary']).read_bytes()).hexdigest()!=state['binary_sha256'] or \
    manifest['binary_sha256']!=state['binary_sha256']:
    error('Controlled binary does not match frozen candidate')
for item in manifest['source_files']:
    if hashlib.sha256((candidate/'source'/item['path']).read_bytes()).hexdigest()!=item['sha256']:
        error('Frozen source snapshot changed: '+item['path'])
for name,expected_hash in manifest['test_evidence_sha256'].items():
    if hashlib.sha256((candidate/'tests'/name).read_bytes()).hexdigest()!=expected_hash:
        error('Frozen test evidence changed: '+name)

product_epochs={}
for marker in sorted((run/'outputs/zhang_internal_products.csv.epochs').glob('*/COMMIT.json')):
    commit=json.loads(marker.read_text())
    if commit.get('ar') is not True:continue
    with (marker.parent/'products.csv').open() as stream:
        products=[r for r in csv.DictReader(stream) if r['solution']=='PRODUCT_FIXED']
    if not products:error('No fixed product rows at '+commit['epoch'])
    product_epochs[commit['epoch'][-8:]]={(r['satellite'],r['observable']):r for r in products}
    if not any(r['pppar_usable']=='1' for r in products):
        error('No usable fixed product at '+commit['epoch'])
    max_mean_error=max_cov_error=0.
    for row in products:
        if row['pppar_usable']!='1':continue
        values={k:float(row[k]) for k in ('clock_m','phase_m','correction_m','clock_sigma_m',
            'phase_sigma_m','correction_sigma_m','clock_phase_covariance_m2')}
        if not all(math.isfinite(v) for v in values.values()):error('Nonfinite usable product '+commit['epoch'])
        mean_error=abs(values['correction_m']-(values['clock_m']-values['phase_m']))
        cov_error=abs(values['correction_sigma_m']**2-(values['clock_sigma_m']**2+
            values['phase_sigma_m']**2-2*values['clock_phase_covariance_m2']))
        max_mean_error=max(max_mean_error,mean_error);max_cov_error=max(max_cov_error,cov_error)
        if mean_error>1e-7 or cov_error>1e-11:error('Joint correction identity failed '+commit['epoch'])
    covariance={}
    with (marker.parent/'covariance.csv').open() as stream:
        for row in csv.DictReader(stream):
            if row['solution']!='PRODUCT_FIXED':continue
            key=(row['row_satellite'],row['row_parameter'],row['row_observable'])
            other=(row['column_satellite'],row['column_parameter'],row['column_observable'])
            value=float(row['covariance_m2'])
            if not math.isfinite(value):error('Nonfinite joint covariance '+commit['epoch'])
            covariance[key,other]=value
    required=set()
    for row in products:
        if row['pppar_usable']=='1':
            required.update(((row['satellite'],'CLOCK','NONE'),(row['satellite'],'PHASE',row['observable'])))
    if any((a,b) not in covariance for a in required for b in required):
        error('Incomplete usable-product joint covariance '+commit['epoch'])
    if any(abs(covariance[a,b]-covariance[b,a])>1e-11 for a in required for b in required
           if (a,b) in covariance and (b,a) in covariance):error('Asymmetric joint covariance '+commit['epoch'])
    report['epochs'].append(dict(epoch=commit['epoch'],products=len(products),usable=sum(r['pppar_usable']=='1' for r in products),
        resets=sum(r['requires_user_ambiguity_reset']=='1' for r in products),
        mean_identity_max_m=max_mean_error,covariance_identity_max_m2=max_cov_error))

expected_epoch=datetime(2024,7,17,2)
end_epoch=datetime.fromisoformat('2024-07-17 '+('02:59:30' if completed_hour else state['stop_at']))
expected_epochs=[]
while expected_epoch<=end_epoch:
    expected_epochs.append(expected_epoch.strftime('%H:%M:%S'))
    expected_epoch+=timedelta(seconds=30)
if set(product_epochs)!=set(expected_epochs) or len(report['epochs'])!=len(expected_epochs):
    error('Controlled replay does not cover every configured 30-second epoch through stop target')

tokens=('ZHANG_SHARED_PRODUCT_EVENT','ZHANG_CURRENT_EVENT','ZHANG_BRIDGE_',
    'ZHANG_ROUTE_CONTINUITY_SUPPORT','ZHANG_PRODUCT_HANDOFF','ZHANG_HANDOFF_FLOAT',
    'ZHANG_PRODUCT_PUBLICATION_GATE','ZHANG_PRODUCT_NUMERIC_GATE')
selected={epoch:[] for epoch in ('02:03:00','02:03:30','02:04:00','02:06:00','02:37:00','02:37:30')}
policy=epoch_only=0
trace=run/'outputs/Network-case-202419900.TRACE'
with trace.open(errors='replace') as stream:
    for line in stream:
        if 'R51_VALIDATION_POLICY' in line:
            policy+=1
            if not re.search(r'ratio_threshold=2(?:\s|$)',line) or 'nis_gate=0' not in line:
                error('Unexpected runtime AR acceptance policy')
        if 'ZHANG_R51_EPOCH_ONLY_AR ' in line:epoch_only+=1
        if any(token in line for token in tokens):
            for epoch in selected:
                if 'time=2024-07-17 '+epoch in line:selected[epoch].append(line.rstrip())
if not policy or not epoch_only:error('Missing runtime ratio/NIS/epoch-only evidence')
(evidence/'event_trace.txt').write_text('\n'.join(line for lines in selected.values() for line in lines))
committed={row['epoch'][-8:] for row in report['epochs']}
for epoch in ('02:06:00','02:37:30'):
    text='\n'.join(selected[epoch])
    forbidden=('SHARED_EVENT_TOTAL_COOPERATIVE_BUDGET_EXHAUSTED',
        'SHARED_EVENT_PREPARATION_COOPERATIVE_BUDGET_EXHAUSTED','SHARED_EVENT_FINAL_QUERY_COOPERATIVE_BUDGET_EXHAUSTED')
    passed=epoch in committed and bool(text) and not any(word in text for word in forbidden)
    previous={'02:06:00':'02:05:30','02:37:30':'02:37:00'}[epoch]
    current=product_epochs.get(epoch,{})
    prior=product_epochs.get(previous,{})
    common=[key for key in prior if prior[key]['pppar_usable']=='1']
    failed_continuity=[key for key in common if key not in current or current[key]['pppar_usable']!='1' or
        current[key]['same_segment_continuous']!='1' or current[key]['requires_user_ambiguity_reset']!='0' or
        current[key]['public_phase_segment_id']!=prior[key]['public_phase_segment_id']]
    # A positive bridge/publication outcome is necessary; absence of an error is not a proof.
    passed=passed and bool(common) and not failed_continuity
    passed=passed and any('route=SAME_SEGMENT' in line for line in selected[epoch])
    if epoch=='02:37:30':
        passed=passed and 'CURRENT_CONSTRAINTS_CANNOT_ELIMINATE_UNAVAILABLE_ARC' not in text
        mechanisms=[name for name,found in (
            ('CORRELATED_OLD_ARC_FLOAT',any('AGGO|G10|V0' in line and 'ZHANG_BRIDGE_AUXILIARY_SUPPORT' in line
                for line in selected[epoch])),
            ('CURRENT_EXACT_REWRITE','ZHANG_CURRENT_EVENT_REPRESENTATION_REWRITE' in text),
            ('VERIFIED_CURRENT_ROUTE',any('ZHANG_ROUTE_CONTINUITY_SUPPORT' in line and
                ('float_valid=1 float_missing_rank=0' in line or
                 'current_network_valid=1 current_network_missing_rank=0' in line) for line in selected[epoch]))) if found]
        passed=passed and bool(mechanisms)
        passed=passed and sum(key[0]=='G10' for key in common)==2
    timeout_lines=[line for line in selected[epoch] if any(word in line for word in
        ('WALL_TIMEOUT','COOPERATIVE_BUDGET_EXHAUSTED','EXACT_OPERATION_LIMIT','SEARCH_NOT_COMPLETED'))]
    report['events'][epoch]=dict(verified=passed,evidence_lines=len(selected[epoch]),
        compared_products=len(common),failed_continuity=failed_continuity,
        intermediate_timeout_or_unfinished_lines=timeout_lines)
    if epoch=='02:37:30':report['events'][epoch]['evidence_mechanisms']=mechanisms
    if not passed:error('Event did not satisfy evidence gate: '+epoch)

# Known G19 false-same-segment regression. This is a conservative excursion
# alarm, NOT an absolute accuracy or real-time performance certification.
report['g19_same_segment_regression']=[]
for before,after in (('02:03:00','02:03:30'),('02:03:30','02:04:00')):
    prior=product_epochs.get(before,{})
    current=product_epochs.get(after,{})
    for key,row in current.items():
        if key[0]!='G19' or key not in prior:continue
        reference=prior[key]
        if row['pppar_usable']!='1' or reference['pppar_usable']!='1':
            error('G19 regression lacks usable comparison '+after+' '+key[1]);continue
        steps=[float(value['correction_m'])-float(prior[other]['correction_m'])
            for other,value in current.items() if other in prior and other[0]!='G19' and other[1]==key[1] and
            value['pppar_usable']=='1' and prior[other]['pppar_usable']=='1' and
            value['same_segment_continuous']=='1' and
            value['public_phase_segment_id']==prior[other]['public_phase_segment_id']]
        if len(steps)<5:
            error('G19 regression lacks five independent same-segment references '+after+' '+key[1]);continue
        residual=float(row['correction_m'])-float(reference['correction_m'])-statistics.median(steps)
        same=row['same_segment_continuous']=='1' and row['public_phase_segment_id']==reference['public_phase_segment_id']
        report['g19_same_segment_regression'].append(dict(epoch=after,observable=key[1],
            common_mode_removed_step_m=residual,same_segment=same,excursion_alarm_m=0.5))
        if same and abs(residual)>0.5:error('G19 unsupported same-segment excursion '+after+' '+key[1])
if len(report['g19_same_segment_regression'])!=4:error('Incomplete G19 known-event regression')

driver=candidate/'zhang_r51_bridge_event_replay_tests'
if hashlib.sha256(driver.read_bytes()).hexdigest()!=manifest['exact_replay_driver_sha256']:
    error('Frozen exact replay driver changed')
inputs=evidence/'bridge_inputs'
if not list(inputs.glob('*.txt')):error('No real pre-image bridge captures')
else:
    with (evidence/'exact_event_replay.jsonl').open('w') as stream:
        result=subprocess.run([str(driver),str(inputs)],stdout=stream,stderr=subprocess.STDOUT)
    if result.returncode:error('Exact old/new real-event replay failed')
    else:
        results=[json.loads(line) for line in (evidence/'exact_event_replay.jsonl').read_text().splitlines() if line.startswith('{')]
        if not results or any(row.get('exact_equal') is not True for row in results):
            error('Missing positive exact old/new image comparisons')
        for epoch in ('02:06:00','02:37:30'):
            if not any(row.get('epoch','').endswith(epoch) for row in results):
                error('No real exact replay coverage at '+epoch)
        if not any(row.get('epoch','').endswith('02:37:30') and row.get('satellite')=='G10' for row in results):
            error('No G10 exact image replay at old-arc event')

unit_xml=candidate/'tests/bridge_ctest.xml'
try:
    root=ET.parse(unit_xml).getroot()
    cases=list(root.iter('testcase'))
    if int(root.attrib.get('failures','0')) or int(root.attrib.get('errors','0')) or not cases or any(
        case.find('failure') is not None or case.find('error') is not None or case.find('skipped') is not None
        for case in cases):
        error('CTest XML contains failures or no tests')
    long_names={'zhang_cycle_slip_'+name for name in ('no_slip_holdout','injected_noise_readiness',
        'doppler_no_slip','doppler_readiness','storm_no_slip','storm_injected','doppler_blunder',
        'smooth_phase_trend','no_doppler_iono_trend','no_doppler_phase_trend')}
    configured=json.loads((candidate/'tests/configured_tests.json').read_text())
    expected={test['name'] for test in configured['tests']}-long_names
    if {case.attrib['name'] for case in cases}!=expected:error('CTest XML does not cover all required non-long tests')
except (OSError,ET.ParseError):error('No completed CTest XML')
long_log=candidate/'tests/cycle_slip_unchanged_regression.log'
if '100% tests passed, 0 tests failed out of 10' not in long_log.read_text():
    error('Long cycle-slip regressions are incomplete or failed')
qc_result=Path('/mnt/c/Users/rx/Documents/GINAN/artifacts/cycle_slip_current_epoch_20261018/validation_v3/validation_result.json')
qc=json.loads(qc_result.read_text())
repo=candidate/'source'
qc_files=('src/cpp/common/acsQC.cpp','src/cpp/common/zhangCycleSlipEvidence.hpp',
    'src/cpp/common/zhangCycleSlipDopplerWindow.hpp','src/cpp/common/zhangCycleSlipEvidenceState.hpp')
if qc.get('all_execution_checks_passed') is not True or any(
    hashlib.sha256((repo/path).read_bytes()).hexdigest()!=qc['source_hashes'][path] for path in qc_files):
    error('Independent QC readiness evidence no longer matches unchanged detector sources')
performance=qc.get('performance_gates',[])
if {item.get('name') for item in performance}!={'injected_independent_holdout','storm_injected_holdout',
    'without_doppler_holdout'} or any(item.get('passed') is not True for item in performance):
    error('Independent cycle-slip performance gates are missing or failed')
report['unchanged_detector_independent_readiness']=dict(source=str(qc_result),
    performance_gates=qc.get('performance_gates'),reuse_requires_exact_source_hashes=True)

checkpoint_driver=candidate/'structural_checkpoint_roundtrip'
if hashlib.sha256(checkpoint_driver.read_bytes()).hexdigest()!=manifest['checkpoint_replay_driver_sha256']:
    error('Frozen structural checkpoint driver changed')
new_checkpoints=sorted((run/'checkpoints').glob('*/checkpoint.bundle'))
if not new_checkpoints:error('No new schema6 checkpoint for real production roundtrip')
else:
    with (evidence/'structural_checkpoint_roundtrip.log').open('w') as stream:
        result=subprocess.run([str(checkpoint_driver),str(new_checkpoints[-1]),state['binary_sha256']],
            stdout=stream,stderr=subprocess.STDOUT)
    if result.returncode:error('Real schema6 production structural checkpoint roundtrip failed')

gates=dict(binary_sha256=state['binary_sha256'],all_required_tests_passed=not report['errors'],
    event_0206_verified=report['events']['02:06:00']['verified'],
    event_0237_verified=report['events']['02:37:30']['verified'],
    ratio=2,joint_nis=False,epoch_only=True,ar_started=completed_hour,
    completed_hour_mode=completed_hour,errors=report['errors'])
for key in ('config_sha256','pinned_environment','checkpoint_sha256','checkpoint_manifest_sha256','launcher_sha256'):
    gates[key]=state[key]
(evidence/(report_name+'.json')).write_text(json.dumps(report,indent=2))
temporary=evidence/(gate_name+'.tmp')
temporary.write_text(json.dumps(gates,indent=2))
temporary.replace(evidence/(gate_name+'.json'))
print(json.dumps(gates,indent=2),flush=True)
sys.exit(0 if gates['all_required_tests_passed'] else 1)
