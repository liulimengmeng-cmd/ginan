"""Replay every real WL/L1 trial through separate frozen/current executables.
WORK_ROOT. Run only after both search stages have closed their input streams.
"""
from pathlib import Path
import json
import subprocess
import sys

work=Path(sys.argv[1]).resolve();repo=Path(__file__).resolve().parents[2]
selected=json.loads((work/'candidate_selection.json').read_text()) if (work/'candidate_selection.json').exists() else {}
suffix='_'+selected['tag'] if selected.get('tag') else ''
probe=work/selected.get('probe_directory','one_epoch_probe')
results={}
for stage in ('WL','L1'):
    source=probe/'condition_inputs'/(stage+'.bin')
    prefix=stage+suffix
    for mode,folder in (('old','baseline_condition'),('new','current')):
        output=work/(prefix+'_'+mode+'.posterior')
        with (work/(prefix+'_'+mode+'_replay.log')).open('w') as log:
            subprocess.run([str(work/folder/'incremental_condition_contract'),'replay',str(source),str(output)],
                stdout=log,stderr=subprocess.STDOUT,check=True)
    result=work/(prefix+'_replay_comparison.json')
    with (work/(prefix+'_replay_comparison.log')).open('w') as log:
        subprocess.run(['python3',str(repo/'src/tests/compare_r51_conditioners.py'),
            str(work/(prefix+'_old.posterior')),str(work/(prefix+'_new.posterior')),str(result)],
            stdout=log,stderr=subprocess.STDOUT,check=True)
    results[stage]=json.loads(result.read_text());print(stage,results[stage],flush=True)
(work/selected.get('real_replay_report','real_conditioner_replays.json')).write_text(json.dumps(results,indent=2))
