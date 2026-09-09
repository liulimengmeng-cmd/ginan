#!/usr/bin/env python3
"""Read completed Experiment 5 runs; write new audit artifacts without replacement."""
import argparse
from collections import Counter
import json
from pathlib import Path
import re
from experiment5_code_osb_sigma import sha, write_json, DOC, C
import audit_external_bias_application as bias
import audit_pppar_feedback_information as feedback
import audit_pppar_dual_frequency_datum as datum
import audit_pppar_coordinate_integrity as coordinate


def network(directory):
    paths = [p for p in Path(directory).glob('Network-*.trace') if '_smoothed' not in p.name]
    assert len(paths) == 1, paths
    return paths[0]


def screening(trace):
    records = []
    current = None
    epoch = None
    with trace.open(errors='replace') as f:
        for line in f:
            m = re.search(r'Epoch\s+(\d+)\s+=',line)
            if m:
                epoch = int(m[1])
            if 'PPP_AR FEEDBACK_SHADOW_SUMMARY ' in line:
                current = {'epoch':epoch, 'postfit_failure_count':0, 'prefit_failure_count':0,
                           'iteration_limit_count':0, 'filter_failed_count':0}
            if current is not None:
                for text, key in [('Postfit check failed','postfit_failure_count'),
                                  ('Prefit check failed','prefit_failure_count'),
                                  ('iterations limit reached','iteration_limit_count'),
                                  ('FILTER FAILED','filter_failed_count')]:
                    if text in line:
                        current[key] += 1
                if 'PPP_AR FEEDBACK_SHADOW_REALISATION ' in line:
                    values = feedback.fields(line)
                    current['mixed_unit_shadow_norm'] = float(values.get('state_update_error_norm','nan'))
                    records.append(current)
                    current = None
    with_screening = [r for r in records if r['postfit_failure_count'] or r['prefit_failure_count']]
    without_screening = [r for r in records if not r['postfit_failure_count'] and not r['prefit_failure_count']]
    return {'feedback_calls': len(records), 'calls_with_screening_failure':len(with_screening),
            'postfit_failure_events':sum(r['postfit_failure_count'] for r in records),
            'filter_failed_events':sum(r['filter_failed_count'] for r in records),
            'iteration_limit_events':sum(r['iteration_limit_count'] for r in records),
            'mixed_unit_norm_with_screening': feedback.summarize([r['mixed_unit_shadow_norm'] for r in with_screening]),
            'mixed_unit_norm_without_screening': feedback.summarize([r['mixed_unit_shadow_norm'] for r in without_screening]),
            'interpretation': 'Full-state norm mixes metres, cycles and TECU; it is not a coordinate error in metres.',
            'largest_discrepancy_calls': sorted(records,key=lambda r:r['mixed_unit_shadow_norm'],reverse=True)[:10]}


def prior_check(directory, expected_sigma):
    inits = Counter()
    variances = set()
    incorrect_mode = per_epoch_nonzero = 0
    for path in Path(directory).glob('*.trace'):
        if path.name.startswith('Network-'):
            continue
        with path.open(errors='replace') as f:
            for line in f:
                if 'PPP_EXTERNAL_BIAS_APPLICATION type=PHASE ' not in line:
                    continue
                m = bias.MARKER_RE.search(line)
                assert m
                d = m.groupdict()
                incorrect_mode += d['variance_mode'] != 'PERSISTENT_STATE'
                per_epoch_nonzero += float(d['applied_variance_m2']) != 0
                variances.add(float(d['state_prior_variance_m2']))
                if d['state_prior_initialised'] == '1':
                    inits[d['satellite']+'/'+d['signal']] += 1
    expected = (float(f'{expected_sigma/C*1e9:.7f}')*C/1e9)**2
    return {'initialisations_by_satellite_signal':dict(sorted(inits.items())),
            'distinct_prior_variances_m2':sorted(variances), 'expected_encoded_prior_variance_m2':expected,
            'incorrect_mode_count':incorrect_mode, 'phase_per_epoch_variance_nonzero_count':per_epoch_nonzero,
            'each_of_62_states_initialised_once':len(inits)==62 and set(inits.values())=={1},
            'prior_variance_matches_fixture':bool(variances) and all(abs(v-expected)<1e-16 for v in variances)}


def main(root):
    plan = json.loads((root/'plan.json').read_text())
    audit = root/'audits'
    audit.mkdir(exist_ok=False)
    for run in plan['runs']:
        label = run['label']
        completed = json.loads((root/f'{label}.finished.json').read_text())
        assert completed['returncode'] == 0
        trace = network(run['output'])
        reports = {'feedback':feedback.audit_trace(trace),
                   'bias':bias.audit_trace(Path(run['output'])),
                   'datum':datum.audit_trace(trace),
                   'prior':prior_check(run['output'],run['assumed_sigma_m']),
                   'screening':screening(trace)}
        for name,value in reports.items():
            write_json(audit/f'{label}_{name}.json',value)
        assert reports['feedback']['epoch_count'] == 2880
        assert reports['bias']['duplicate_application_key_count'] == 0
        assert reports['bias']['missing_external_bias_count'] == 0
        assert reports['prior']['each_of_62_states_initialised_once']
        assert reports['prior']['prior_variance_matches_fixture']
        if run['mode']=='float':
            assert reports['feedback']['feedback_epoch_count']==0
        print('audited',label,flush=True)
    for mm in (3,10):
        result = coordinate.audit_coordinate_integrity(
            network(root/f'sigma_{mm}mm_ar'),network(root/f'sigma_{mm}mm_float'),
            Path('/home/rx/GINAN/inputData/products/IGS0OPSSNX_20241960000_07D_07D_CRD.SNX'),
            coordinate.Thresholds(), primary_state_block='AR',float_state_block='PPP',
            comparison_burn_in_minutes=840)
        write_json(audit/f'sigma_{mm}mm_coordinate.json',coordinate.compact_report(result))
    write_json(audit/'baseline_zero_screening.json',screening(network('/mnt/d/tec/code-exp4c-full-20260909')))
    originals = json.loads((root/'all_inputs.json').read_text())
    for f in originals:
        assert sha(f['absolute_path'])==f['sha256']
    write_json(audit/'original_input_preservation.json',{'verified_sha256_count':len(originals),'all_unchanged':True})


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('root',type=Path)
    main(p.parse_args().root)
