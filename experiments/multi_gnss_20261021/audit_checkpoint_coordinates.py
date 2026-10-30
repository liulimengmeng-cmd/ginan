#!/usr/bin/env python3
"""Read-only static XYZ/SNX consistency, NOT independent accuracy validation.

Input inventory comes from the authenticated C++ checkpoint reader. Expected
stations come from the run's pinned observations, not from available estimates.
"""
import argparse
from datetime import datetime, timedelta
import hashlib
import json
import math
from pathlib import Path
import statistics

AXES = ('X', 'Y', 'Z')
EXPECTED_MODE = 'E29_MULTI_GNSS_FLOAT_V2_GPS_GAL_BDS'

def snx_epoch(value):
    year, day, second = map(int, value.split(':'))
    if year == day == second == 0:
        return None
    year += 1900 if year >= 80 else 2000
    return datetime(year, 1, 1) + timedelta(days=day - 1, seconds=second)

def sinex_positions(path):
    groups, intervals = {}, {}
    block = ''
    for line in path.read_text(errors='strict').splitlines():
        if line.startswith('+'):
            block = line[1:].strip()
        elif line.startswith('-'):
            block = ''
        elif line.startswith(' ') and not line.startswith('*'):
            fields = line.split()
            if block == 'SOLUTION/ESTIMATE' and len(fields) >= 10 and fields[1] in (
                    'STAX', 'STAY', 'STAZ', 'VELX', 'VELY', 'VELZ'):
                key = tuple(fields[2:5])
                item = {'epoch': fields[5], 'unit': fields[6], 'value': float(fields[8].replace('D', 'E'))}
                values = groups.setdefault(key, {})
                if fields[1] in values:
                    raise ValueError('Duplicate SINEX solution parameter: ' + str(key))
                values[fields[1]] = item
            elif block == 'SOLUTION/EPOCHS' and len(fields) >= 7:
                key = tuple(fields[:3])
                if key in intervals:
                    raise ValueError('Duplicate SINEX solution interval: ' + str(key))
                intervals[key] = (snx_epoch(fields[4]), snx_epoch(fields[5]))
    return groups, intervals

def reference_for(station, epoch, groups, intervals):
    candidates = []
    for key, values in groups.items():
        if key[0] != station or not all('STA' + axis in values for axis in AXES):
            continue
        start, end = intervals.get(key, (None, None))
        if start and epoch < start or end and epoch > end:
            continue
        candidates.append((key, values))
    if not candidates:
        return {'status': 'NO_REFERENCE_SOLUTION'}
    if len(candidates) != 1:
        return {'status': 'AMBIGUOUS_REFERENCE_SOLUTION', 'candidates': [list(k) for k, _ in candidates]}
    key, values = candidates[0]
    epochs = {values['STA' + axis]['epoch'] for axis in AXES}
    if len(epochs) != 1 or any(values['STA' + axis]['unit'] != 'm' for axis in AXES):
        return {'status': 'INCONSISTENT_REFERENCE_EPOCH_OR_UNITS'}
    ref_text = epochs.pop()
    ref_epoch = snx_epoch(ref_text)
    if ref_epoch is None:
        return {'status': 'UNDEFINED_REFERENCE_EPOCH'}
    seconds = (epoch - ref_epoch).total_seconds()
    xyz = [values['STA' + axis]['value'] for axis in AXES]
    velocity = all('VEL' + axis in values and values['VEL' + axis]['unit'] == 'm/y' for axis in AXES)
    if velocity:
        if any(values['VEL' + axis]['epoch'] != ref_text for axis in AXES):
            return {'status': 'POSITION_VELOCITY_EPOCH_MISMATCH'}
        xyz = [x + values['VEL' + axis]['value'] * seconds / (365.25 * 86400)
               for axis, x in zip(AXES, xyz)]
    return {'status': 'AVAILABLE', 'xyz_m': xyz, 'point': key[1], 'solution': key[2],
            'reference_epoch': ref_text, 'epoch_difference_seconds': seconds,
            'velocity_propagated': velocity,
            'epoch_scope': 'PROPAGATED_USING_SINEX_VELOCITY' if velocity else
                           'UNPROPAGATED_REFERENCE_EPOCH_COORDINATE_NOT_EXACT_TARGET_EPOCH'}

def enu(delta, reference):
    x, y, z = reference
    longitude = math.atan2(y, x)
    p = math.hypot(x, y)
    e2 = 6.6943799901413165e-3
    latitude = math.atan2(z, p * (1 - e2))
    for _ in range(10):
        n = 6378137 / math.sqrt(1 - e2 * math.sin(latitude) ** 2)
        latitude = math.atan2(z + e2 * n * math.sin(latitude), p)
    sl, cl = math.sin(longitude), math.cos(longitude)
    sp, cp = math.sin(latitude), math.cos(latitude)
    dx, dy, dz = delta
    return [-sl * dx + cl * dy, -sp * cl * dx - sp * sl * dy + cp * dz,
            cp * cl * dx + cp * sl * dy + sp * dz]

def summary(values):
    if not values:
        return {'count': 0, 'status': 'UNAVAILABLE'}
    ordered = sorted(abs(x) for x in values)
    return {'count': len(values), 'rms_m': math.sqrt(sum(x*x for x in values) / len(values)),
            'median_m': statistics.median(values), 'p95_absolute_m': ordered[math.ceil(.95 * len(ordered)) - 1],
            'maximum_absolute_m': ordered[-1]}

def audit(inventory, plan, snx):
    if inventory['producer_sha256'] != plan['binary_sha256']:
        raise ValueError('Inventory and pinned run producer hashes differ')
    if inventory.get('experiment_mode') != EXPECTED_MODE:
        raise ValueError('Not an authenticated target G/E/C FLOAT checkpoint mode')
    expected = [row['station'] for row in plan['pinned_input_contract']['observations']]
    if len(set(expected)) != len(expected):
        raise ValueError('Duplicate expected station ID')
    epoch = datetime.fromisoformat(inventory['checkpoint_epoch'])
    groups, intervals = sinex_positions(snx)
    available = inventory['coordinates']
    results, e, n, u, three = {}, [], [], [], []
    for station in expected:
        row = {'station': station}
        estimate = available.get(station)
        if estimate is None:
            row['status'] = 'NO_REC_POS_STATE'
        elif not estimate['complete_xyz'] or not estimate['finite_xyz_positive_variance'] or not estimate.get('absolute_ecef_plausible', False):
            row.update(status='INCOMPLETE_OR_INVALID_REC_POS', estimate=estimate)
        else:
            reference = reference_for(station, epoch, groups, intervals)
            row.update(estimate=estimate, reference=reference, status=reference['status'])
            if reference['status'] == 'AVAILABLE':
                delta = [x - r for x, r in zip(estimate['xyz_m'], reference['xyz_m'])]
                east, north, up = enu(delta, reference['xyz_m'])
                distance = math.sqrt(sum(x*x for x in delta))
                row.update(enu_difference_m=[east, north, up], xyz_difference_m=delta, difference_3d_m=distance)
                e.append(east); n.append(north); u.append(up); three.append(distance)
        results[station] = row
    complete = all(available.get(s, {}).get('complete_xyz', False) and
                   available.get(s, {}).get('finite_xyz_positive_variance', False) and
                   available.get(s, {}).get('absolute_ecef_plausible', False) and
                   available.get(s, {}).get('coordinate_block_positive_definite', False) for s in expected)
    return {'schema': 'STATIC_FLOAT_COORDINATE_CONSISTENCY_V1', 'checkpoint_epoch': inventory['checkpoint_epoch'],
            'producer_sha256': inventory['producer_sha256'], 'bundle_sha256': inventory['bundle_sha256'],
            'is_planned_final_epoch': inventory['checkpoint_epoch'] == plan.get('end_epoch'),
            'experiment_mode': inventory['experiment_mode'],
            'planned_epoch_count': plan.get('expected_epoch_count'),
            'span_evidence_scope': 'CHECKPOINT_ALONE_IS_NOT_PROOF_OF_ALL_241_COMMITTED_EPOCHS_REQUIRES_RUN_LOG_AUDIT',
            'expected_station_count': len(expected), 'estimated_station_count': len(available),
            'all_expected_XYZ_positive_finite_coordinate_covariance': complete,
            'unexpected_estimated_stations': sorted(set(available) - set(expected)),
            'reference_available_count': len(three), 'station_results': results,
            'distribution': dict(E=summary(e), N=summary(n), U=summary(u), distance_3d=summary(three)),
            'full_x_finite': inventory['full_x_finite'], 'full_P_finite': inventory['full_P_finite'],
            'full_P_validation_scope': inventory.get('P_validation_scope', 'UNAVAILABLE'),
            'full_P_maximum_asymmetry_mixed_native_units': inventory.get('full_P_maximum_asymmetry_mixed_native_units'),
            'full_P_minimum_diagonal_mixed_native_units': inventory.get('full_P_minimum_diagonal_mixed_native_units'),
            'receiver_position_rate_count': inventory['receiver_position_rate_count'],
            'state_counts_by_system_and_type': inventory['state_counts_by_system_and_type'],
            'reference_sha256': hashlib.sha256(snx.read_bytes()).hexdigest(), 'reference_path': str(snx),
            'reference_claim': 'SAME_SOURCE_WEEKLY_SNX_INITIAL_PRIOR_CONSISTENCY_NOT_INDEPENDENT_ACCURACY',
            'initial_prior_ecef_sigma_m': [1, 1, 1], 'static_coordinate_process_noise': [0, 0, 0],
            'no_helmert_alignment_or_epoch_offsets_fit': True, 'scientific_accuracy_accepted': False}

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--inventory', type=Path, required=True)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--sinex', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = audit(json.loads(args.inventory.read_text()), json.loads(args.plan.read_text()), args.sinex)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({key: result[key] for key in ('expected_station_count', 'estimated_station_count',
          'reference_available_count', 'distribution', 'reference_claim')}, indent=2))
