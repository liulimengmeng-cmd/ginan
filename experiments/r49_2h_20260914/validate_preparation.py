"""Validate completion gates and unchanged scientific settings without PEA."""
from pathlib import Path
import datetime, hashlib, importlib.util, json, re, tempfile
import yaml
W=Path(__file__).resolve().parent
module=importlib.util.spec_from_file_location('gate',W/'check_and_launch.py')
gate=importlib.util.module_from_spec(module);module.loader.exec_module(gate)
spec=json.loads((W/'experiment.json').read_text())
parent=Path(spec['continuation']['parent_audit'])
pspec=json.loads((parent/'experiment.json').read_text())
start=int(datetime.datetime(2024,7,17,tzinfo=datetime.timezone.utc).timestamp())
with tempfile.TemporaryDirectory() as directory:
    fixture=Path(directory);state=fixture/(pspec['case']+'.state.txt')
    report=fixture/'r49_final_verification.json'
    assert gate.completion(fixture,pspec)=='PENDING'
    state.write_text('RUNNER_EXIT status=109\n')
    assert gate.completion(fixture,pspec)=='FAILED'
    state.write_text('RUNNER_EXIT status=0\n')
    assert gate.completion(fixture,pspec)=='FAILED'
    value={'case':pspec['case'],'pea_status':0,'runtime_pass':True,'failures':[],'alarms':{},'epochs':[{'epoch':start+30*i} for i in range(61)]}
    report.write_text(json.dumps(value));assert gate.completion(fixture,pspec)=='PASSED'
    value['alarms']={'DGEMV_ILLEGAL_VALUE':1}
    report.write_text(json.dumps(value));assert gate.completion(fixture,pspec)=='FAILED'
    value['alarms']={};value['epochs'][-1]['epoch']-=30
    report.write_text(json.dumps(value));assert gate.completion(fixture,pspec)=='FAILED'
old=yaml.safe_load((parent/'config/case.yaml').read_text())
new=yaml.safe_load((W/'config/case.yaml').read_text())
for config in [old,new]:
    config['processing_options']['epoch_control'].pop('end_epoch')
    z=config['processing_options']['gnss_general']['zhang_pppar']
    for key in ['checkpoint_output_directory','checkpoint_capture_epochs','product_filename','product_covariance_filename']:z.pop(key)
    config.pop('outputs')
assert old==new,'Scientific settings changed'
inputs=yaml.safe_load((W/'config/zhang_global_2024199_180_inputs.yaml').read_text())
rinex=inputs['inputs']['gnss_observations']['rnx_inputs']
coverage=[]
for name in rinex:
    p=Path('/mnt/d/GINAN_R20/inputData/data')/name
    with p.open('rb') as stream:
        header=stream.read(65536).decode(errors='replace')
        stream.seek(max(0,p.stat().st_size-262144));tail=stream.read().decode(errors='replace')
    matches=re.findall(r'^>\s+(\d{4})\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)',tail,re.M)
    assert matches,'No final RINEX epoch: '+name
    last=datetime.datetime(*map(int,matches[-1]))
    assert last>=datetime.datetime(2024,7,17,2),name
    coverage.append({'station':name[:4],'last_epoch':str(last)})
assert len(coverage)==180
result={'gate_cases_passed':6,'scientific_settings_unchanged':True,'stations_covering_end':len(coverage),'input_tail_coverage':coverage,'note':'Tail coverage does not establish gap-free observations; full input hashes are checked by the runner.'}
(W/'preparation_validation.json').write_text(json.dumps(result,indent=2)+'\n')
files=sorted(p for p in W.rglob('*') if p.is_file() and (p.suffix in ['.py','.sh','.yaml'] or p.name in ['experiment.json','frozen_r49.json','INPUTS_SHA256SUMS','preparation_validation.json']))
(W/'PREPARED_SHA256SUMS').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(W))+'\n' for p in files))
print(json.dumps({k:v for k,v in result.items() if k!='input_tail_coverage'},indent=2))
