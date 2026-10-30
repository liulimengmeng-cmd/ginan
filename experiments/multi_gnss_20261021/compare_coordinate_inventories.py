"""Read-only numerical comparison of two authenticated coordinate inventories."""
import argparse
import hashlib
import json
import math
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('first', type=Path)
parser.add_argument('second', type=Path)
parser.add_argument('output', type=Path)
args = parser.parse_args()
a, b = (json.loads(p.read_text()) for p in (args.first, args.second))
for field in ('producer_sha256', 'experiment_mode', 'checkpoint_epoch', 'runtime_id'):
    if a[field] != b[field]:
        raise ValueError('Unmatched comparison identity: ' + field)
if set(a['coordinates']) != set(b['coordinates']):
    raise ValueError('Different coordinate station sets')
rows = {}
for station, left in a['coordinates'].items():
    right = b['coordinates'][station]
    dx = [y-x for x, y in zip(left['xyz_m'], right['xyz_m'])]
    dq = [y-x for x, y in zip(left['covariance_xyz_row_major_m2'], right['covariance_xyz_row_major_m2'])]
    rows[station] = {
        'B_minus_A_xyz_m': dx,
        'delta_position_3d_m': math.sqrt(sum(x*x for x in dx)),
        'B_minus_A_coordinate_covariance_row_major_m2': dq,
        'maximum_absolute_coordinate_covariance_difference_m2': max(map(abs, dq)),
        'coordinate_covariance_difference_frobenius_m2': math.sqrt(sum(x*x for x in dq)),
        'B_minus_A_coordinate_variance_xyz_m2': [dq[i] for i in (0, 4, 8)],
        'A_sigma_xyz_m': left['sigma_xyz_m'], 'B_sigma_xyz_m': right['sigma_xyz_m'],
        'A_coordinate_block_positive_definite': left['coordinate_block_positive_definite'],
        'B_coordinate_block_positive_definite': right['coordinate_block_positive_definite'],
    }
ranked = sorted(rows, key=lambda s: rows[s]['delta_position_3d_m'], reverse=True)
report = {
    'schema': 'ACTUAL_COLD_A_B_COORDINATE_COMPARISON_V1',
    'producer_sha256': a['producer_sha256'], 'checkpoint_epoch': a['checkpoint_epoch'],
    'A_bundle_sha256': a['bundle_sha256'], 'B_bundle_sha256': b['bundle_sha256'],
    'A_inventory_sha256': hashlib.sha256(args.first.read_bytes()).hexdigest(),
    'B_inventory_sha256': hashlib.sha256(args.second.read_bytes()).hexdigest(),
    'station_count': len(rows), 'station_differences': rows,
    'station_order_by_position_difference': ranked,
    'maximum_position_difference_m': max(row['delta_position_3d_m'] for row in rows.values()),
    'minimum_position_difference_m': min(row['delta_position_3d_m'] for row in rows.values()),
    'rms_position_difference_3d_m': math.sqrt(sum(row['delta_position_3d_m']**2 for row in rows.values()) / len(rows)),
    'maximum_coordinate_covariance_entry_difference_m2': max(row['maximum_absolute_coordinate_covariance_difference_m2'] for row in rows.values()),
    'coordinate_blocks_bitwise_numerically_equal': all(row['delta_position_3d_m'] == 0 and row['coordinate_covariance_difference_frobenius_m2'] == 0 for row in rows.values()),
    'scope': 'TWO_INDEPENDENT_12_STATION_COLD_SHORT_REPLAYS_SAME_EPOCH_NOT_A_RESTORE_CAUSAL_TEST_OR_180_STATION_QUALITY',
}
args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
print(json.dumps({key: value for key, value in report.items() if key != 'station_differences'}, indent=2))
for station in ranked:
    row = rows[station]
    print(station, 'delta_xyz_m=', row['B_minus_A_xyz_m'], 'norm_m=', row['delta_position_3d_m'],
          'max_delta_Qxyz_m2=', row['maximum_absolute_coordinate_covariance_difference_m2'])
