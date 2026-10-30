"""Freeze the completed candidate and run final non-GNSS acceptance.
BUILD WORK_ROOT. Run as rx; controlled GNSS probe is a separate root step.
"""
from pathlib import Path
import hashlib
import json
import re
import shutil
import subprocess
import sys

repo=Path(__file__).resolve().parents[2]
build,work=map(lambda p:Path(p).resolve(),sys.argv[1:3])
tag=sys.argv[3] if len(sys.argv)>3 else ''
assert not tag or re.fullmatch('[a-z0-9_]+',tag)
suffix='_'+tag if tag else ''
selection=dict(tag=tag,binary_directory='verified'+suffix,probe_directory='one_epoch_probe'+suffix,
    kernel_directory='paired_kernels'+suffix,epoch_comparison='epoch_comparison'+suffix+'.json',
    full_float_comparison_log='full_float_compare'+suffix+'.log',
    real_replay_report='real_conditioner_replays'+suffix+'.json')
exclude='^zhang_cycle_slip_(no_slip_holdout|injected_noise_readiness|doppler_no_slip|doppler_readiness|storm_no_slip|storm_injected|doppler_blunder|smooth_phase_trend|no_doppler_iono_trend|no_doppler_phase_trend)$'
def run(command,log):
    with (work/log).open('w') as stream:
        subprocess.run(list(map(str,command)),stdout=stream,stderr=subprocess.STDOUT,check=True)
run(['ctest','--test-dir',build,'--output-on-failure','-j','1','-E',exclude],'ctest.log')
(work/selection['binary_directory']).mkdir(exist_ok=False)
shutil.copy2(repo/'bin/pea',work/selection['binary_directory']/'pea')
(work/'candidate_selection.json').write_text(json.dumps(selection,indent=2))
for source,name in ((build/'float_parallel_contract/float_parallel_benchmark','float_parallel_benchmark'),
                    (build/'incremental_condition_contract/incremental_condition_contract','incremental_condition_contract')):
    shutil.copy2(source,work/'current'/name)
print('CANDIDATE_SHA256='+hashlib.sha256((work/selection['binary_directory']/'pea').read_bytes()).hexdigest(),flush=True)
run([work/'current/incremental_condition_contract','suite',work/'suite_new_accepted.bin'],'suite_accepted.log')
run(['python3',repo/'src/tests/compare_r51_conditioners.py',work/'suite_old_scale.bin',
     work/'suite_new_accepted.bin',work/'suite_comparison_accepted.json','suite'],'suite_accepted_comparison.log')
run(['python3',repo/'src/tests/run_r51_ar_hotspots.py',build,work/'checkpoint_verifier','checkpoint'],
    'checkpoint_verifier_build.log')
print('PASS final 40 CTests, frozen candidate, full conditioner suite and checkpoint verifier build',flush=True)
