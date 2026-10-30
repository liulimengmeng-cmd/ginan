"""Verify complete current/reference FLOAT checkpoints and exact inputs.
WORK_ROOT EVIDENCE_DIRECTORY. Does not launch a solver or alter checkpoints.
"""
from pathlib import Path
import hashlib
import json
import re
import subprocess
import sys

work,evidence=map(lambda p:Path(p).resolve(),sys.argv[1:3])
old_evidence=evidence.parent/'product_closure_optimization_20261019'
selected=json.loads((work/'candidate_selection.json').read_text()) if (work/'candidate_selection.json').exists() else {}
probe=work/selected.get('probe_directory','one_epoch_probe')
new=json.loads((probe/'probe.json').read_text())
old=json.loads((old_evidence/'one_epoch_probe/real_epoch_audit.json').read_text())['probe']
assert not new['solver_still_running'] and new['returncode']==130
command=[str(work/'checkpoint_verifier/ar_checkpoint_verify'),new['checkpoint_bundles'][0],
    new['binary_sha256'],old['checkpoint_bundles'][0],old['binary_sha256']]
log_path=work/selected.get('full_float_comparison_log','full_float_compare.log')
with log_path.open('w') as log:
    subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
full_text=log_path.read_text()
line=next(line for line in full_text.splitlines() if line.startswith('REAL_FLOAT_FULL_COMPARE '))
metrics={k:float(v) for k,v in re.findall(r'(\w+)=([0-9.eE+-]+)',line)}
assert metrics['mean_relative_norm']<=2e-10 and metrics['mean_component_scaled']<=2e-10
assert metrics['covariance_relative_frobenius']<=2e-11 and metrics['covariance_entry_scaled']<=2e-11
suffix='_'+selected['tag'] if selected.get('tag') else ''
(work/('full_checkpoint_comparison'+suffix+'.json')).write_text(json.dumps(dict(
    passed=True,new_binary_sha256=new['binary_sha256'],old_binary_sha256=old['binary_sha256'],
    new_checkpoint=new['checkpoint_bundles'][0],old_checkpoint=old['checkpoint_bundles'][0],
    metrics=metrics),indent=2))
def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()
frames=[]
for source in sorted((probe/'exact_inputs').glob('frame_*.txt')):
    reference=old_evidence/'exact_inputs'/source.name
    assert digest(source)==digest(reference),'exact integer system changed: '+source.name
    frames.append(dict(name=source.name,sha256=digest(source)))
assert len(frames)==8
(work/'exact_input_comparison.json').write_text(json.dumps(dict(passed=True,frames=frames),indent=2))
print(log_path.read_text(),flush=True)
print('PASS 8 exact physical-integer constraint/target frames byte-identical',flush=True)
