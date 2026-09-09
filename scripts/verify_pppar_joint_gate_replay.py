#!/usr/bin/env python3
"""Verify pre-update gate decisions and rejected STEC output against snapshots."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import numpy as np
from scipy.stats import chi2
from scipy.special import erfc
from audit_pppar_same_epoch_joint_nis import read_snapshot
from audit_pppar_feedback_information import fields


def verify(root, baseline, mode):
    trace = next((root / f'replay_{mode}').glob('Network-*.trace'))
    events = {}
    epoch = None
    for line in trace.open(errors='replace'):
        if line.startswith('fixAndHoldAmbiguities:'):
            epoch = line.split(': ', 1)[1].strip()[:19]
            events[epoch] = {'submissions': 0}
        if epoch is None:
            continue
        if 'PPP_AR JOINT_GATE ' in line:
            assert 'gate' not in events[epoch]
            events[epoch]['gate'] = fields(line)
        if 'PPP_AR PSEUDOOBS_SUBMISSION ' in line:
            events[epoch]['submissions'] += 1
    gates = {t: e for t, e in events.items() if 'gate' in e}
    passed = []
    rejected = []
    for t, e in gates.items():
        g = e['gate']
        threshold = chi2.isf(erfc(float(g['sigma_threshold']) / np.sqrt(2)), int(g['rows']))
        assert abs(float(g['threshold']) - threshold) < 1e-8
        is_pass = g['status'] == 'PASSED_NOMINAL_SCREEN'
        assert is_pass == (float(g['nis']) <= threshold)
        assert e['submissions'] == int(is_pass), (t, e)
        (passed if is_pass else rejected).append(t)
    snapshots = {}
    for p in sorted((root / f'snapshots_{mode}').glob('*.snapshot')):
        t, m, keys = read_snapshot(p)
        _, old, oldkeys = read_snapshot(baseline / f'snapshots_{mode}' / p.name)
        assert keys == oldkeys
        errors = {name: float(np.max(np.abs(m[name] - old[name]))) for name in m}
        for name in ('H', 'R', 'rhs'):
            assert errors[name] == 0, (t, name, errors[name])
        assert errors['v'] < 1e-6 and errors['P'] < 1e-6 and errors['x'] < 1e-5
        S = m['H'] @ m['P'] @ m['H'].T + m['R']
        nis = float(m['v'].ravel() @ np.linalg.solve(S, m['v'].ravel()))
        assert abs(nis - float(gates[t]['gate']['nis'])) < 1e-6
        cpp = subprocess.check_output([str(root / 'test_joint_gate'), str(p)], text=True).split()
        assert cpp[0] == gates[t]['gate']['status']
        assert abs(float(cpp[1]) - nis) < 1e-6
        snapshots[t] = {'matrices': m, 'keys': keys, 'report': {
            'gate': gates[t]['gate'], 'baseline_prior_max_abs_errors': errors,
            'independent_nis': nis, 'cpp_offline_nis': float(cpp[1]),
            'stec_state_count': 0, 'stec_state_error_max': 0.,
            'stec_variance_error_max': 0., 'stec_covariance_error_max': 0.}}
    cov = next((root / f'replay_{mode}' / 'ionstec').glob('*.STEC.COV'))
    by_tow = {}
    for t, s in snapshots.items():
        hh, mm, ss = map(int, t[11:].split(':'))
        by_tow[259200 + hh * 3600 + mm * 60 + ss] = s
    current = None
    local = {}
    metadata_count = 0
    for row in csv.reader(cov.open()):
        if not row or row[0].startswith('#'):
            continue
        if row[0] == 'META':
            metadata_count += 1
            current = by_tow.get(int(float(row[2])))
            local = {}
            if current:
                current['report']['metadata'] = row
                if current['report']['gate']['status'] != 'PASSED_NOMINAL_SCREEN':
                    assert row[15] == '0' and row[19] == 'JOINT_CONSISTENCY_GATE_REJECTED'
        if current is None or current['report']['gate']['status'] == 'PASSED_NOMINAL_SCREEN':
            continue
        m, keys, report = current['matrices'], current['keys'], current['report']
        if row[0] == 'STATE':
            index = int(row[7]); local[int(row[3])] = index
            assert keys[index]['receiver'] == row[4] and keys[index]['satellite'] == row[5]
            report['stec_state_count'] += 1
            report['stec_state_error_max'] = max(report['stec_state_error_max'], abs(float(row[8])-m['x'][index, 0]))
            report['stec_variance_error_max'] = max(report['stec_variance_error_max'], abs(float(row[9])-m['P'][index, index]))
        if row[0] == 'COV':
            a, b = local[int(row[3])], local[int(row[4])]
            report['stec_covariance_error_max'] = max(report['stec_covariance_error_max'], abs(float(row[5])-m['P'][a, b]))
    assert metadata_count == 2504
    for s in snapshots.values():
        report = s['report']
        if report['gate']['status'] != 'PASSED_NOMINAL_SCREEN':
            assert report['stec_state_count'] > 0
            assert report['stec_state_error_max'] < 1e-12
            assert report['stec_covariance_error_max'] < 1e-12
            assert report['stec_variance_error_max'] < 1e-12
    return {'mode': mode, 'epochs': metadata_count, 'candidate_epochs': len(gates),
            'passed_nominal_screen': len(passed), 'rejected': len(rejected),
            'first_pass_GPST': min(passed) if passed else None,
            'first_reject_GPST': min(rejected) if rejected else None,
            'passed_epochs_GPST': passed, 'rejected_epochs_GPST': rejected,
            'snapshots': {t: s['report'] for t, s in snapshots.items()},
            'interpretation': 'Nominal screen pass is not verified integer correctness; false rejection rate unknown'}


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('root', type=Path)
    p.add_argument('baseline', type=Path)
    p.add_argument('mode', choices=['10mm', '0mm'])
    args = p.parse_args()
    result = verify(args.root, args.baseline, args.mode)
    with (args.root / f'{args.mode}_gate_verification.json').open('x') as f:
        json.dump(result, f, indent=2)
    print(json.dumps({k: v for k, v in result.items() if k not in ('snapshots', 'passed_epochs_GPST', 'rejected_epochs_GPST')}, indent=2))
