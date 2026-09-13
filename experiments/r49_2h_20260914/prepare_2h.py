"""Prepare an independent two-hour R49 run without launching a solver."""
from pathlib import Path
import datetime, hashlib, json, re, shutil
import yaml

W = Path(__file__).resolve().parent
OLD = W.parent / 'r49_work_20260913'
FROZEN = Path('/home/rx/GINAN/frozen-r49-condition-domain-20260913')
spec = json.loads((OLD / 'experiment.json').read_text())
old_case = spec['case']
case = old_case.replace('003000_20260913', '020000_20260914')
checkpoint = '/mnt/d/GINAN_R20/r49_2h_checkpoints_20260914'
snapshots = '/mnt/d/GINAN_R20/r49_2h_root_snapshots_20260914'
spec['case'] = case
spec['observations'].update(end='2024-07-17 02:00:00', expected_epochs=241)
spec['checkpoints'].update(directory=checkpoint, expected_bundles=241)
spec['r49_policy']['root_snapshot_directory'] = snapshots
spec['continuation'] = {'parent_audit': str(OLD), 'mode': 'independent fresh replay',
    'gate': 'parent runner exit 0 and final runtime_pass with 61 epochs and no alarms',
    'unchanged_server_binary': True, 'input_data_copied': False}
spec['claims'] += ' Two-hour extension; matched comparisons use common epochs only.'
(W / 'experiment.json').write_text(json.dumps(spec, indent=2, ensure_ascii=False)+'\n')
(W / 'config').mkdir(exist_ok=False)
for p in (FROZEN / 'config').glob('*.yaml'):
    value = p.read_text().replace(old_case, case).replace(
        '/mnt/d/GINAN_R20/r49_checkpoints_20260913', checkpoint)
    if p.name == 'case.yaml':
        value = value.replace('end_epoch: 2024-07-17 00:30:00', 'end_epoch: 2024-07-17 02:00:00')
        start = datetime.datetime(2024, 7, 17)
        epochs = ''.join('        - "'+str(start+datetime.timedelta(seconds=30*i))+'"\n' for i in range(241))
        value = re.sub(r'(?m)(      checkpoint_capture_epochs:\n)(?:        - .*\n)+', lambda m:m[1]+epochs, value)
    (W / 'config' / p.name).write_text(value)
for name in ['INPUTS_SHA256SUMS', 'metrics.py', 'frozen_r49.json']:
    shutil.copyfile(OLD / name, W / name)
verify = (OLD / 'verify_r49.py').read_text().replace('1405211401','1405216801').replace("!=61", "!=241").replace('/61.', '/241.')
(W / 'verify_r49.py').write_text(verify)
runner = (OLD / 'run_r49_30m.sh').read_text().replace(str(OLD),str(W)).replace(old_case,case)
runner = runner.replace('r49_root_snapshots_20260913','r49_2h_root_snapshots_20260914').replace('r49_checkpoints_20260913','r49_2h_checkpoints_20260914')
runner = runner.replace('expected_epochs=61','expected_epochs=241').replace('/00:30:00','/02:00:00')
runner = runner.replace('"$frozen/config/', '"$audit/config/')
# Capture failures during preflight as well as solver/verification failures.
record_start = runner.index('record()')
record_end = runner.index('python3 "$audit/verify_r49.py"', record_start)
record = runner[record_start:record_end]
runner = runner[:record_start]+runner[record_end:]
anchor = 'exec 8>/home/rx/GINAN/.ginan_global_pea_experiment.lock'
runner = runner.replace(anchor, record+anchor)
# Parent gate and exclusive launch claim are owned by check_and_launch.py.
runner = runner.replace(' && ! -e "$log" && ! -e "$state"', '')
runner = runner.replace('cd "$frozen"', 'cd "$audit"\nsha256sum -c PREPARED_SHA256SUMS >"$audit/prepared_preflight_hashes.log" || exit 93\ncd "$frozen"', 1)
(W / 'run_r49_2h.sh').write_text(runner)
merged = {}
def merge(a,b):
    for k,v in b.items():
        if isinstance(v,dict) and isinstance(a.get(k),dict): merge(a[k],v)
        else: a[k]=v
order=['zhang_global_2024199_180_base.yaml','zhang_global_2024199_180_inputs.yaml','zhang_global_2024199_180_product.yaml','zhang_global_2024199_180_e29_180network_product_1h.yaml','zhang_global_2024199_180_e29_service_30s_1h.yaml','case.yaml']
for name in order: merge(merged,yaml.safe_load((W/'config'/name).read_text()))
control=merged['processing_options']['epoch_control']
assert str(control['end_epoch']) == '2024-07-17 02:00:00' and control['epoch_interval']==30
assert len(merged['processing_options']['gnss_general']['zhang_pppar']['checkpoint_capture_epochs'])==241
assert hashlib.sha256((FROZEN/'bin/pea').read_bytes()).hexdigest()==json.loads((W/'frozen_r49.json').read_text())['binary_sha256']
print(json.dumps({'prepared':str(W),'case':case,'epochs':241,'solver_launched':False},indent=2))
