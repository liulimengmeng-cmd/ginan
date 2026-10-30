import json
import math
from pathlib import Path
import tempfile
import unittest
from datetime import datetime

import audit_checkpoint_coordinates as audit


class CoordinateAuditTests(unittest.TestCase):
    def test_enu_at_equator(self):
        self.assertEqual(audit.enu([1, 2, 3], [6378137, 0, 0]), [2, 3, 1])

    def test_enu_rotation_preserves_norm(self):
        delta = [0.2, -0.5, 0.7]
        result = audit.enu(delta, [1113194, -4841695, 3985350])
        self.assertAlmostEqual(sum(x*x for x in delta), sum(x*x for x in result), places=14)

    def test_zero_and_gpst_epoch(self):
        self.assertIsNone(audit.snx_epoch('00:000:00000'))
        self.assertEqual(audit.snx_epoch('24:199:43200'), datetime(2024, 7, 17, 12))

    def fixture(self, velocity=False):
        rows = ['+SOLUTION/APRIORI',
                ' 1 STAX ABCD A 1 24:199:43200 m 2 9999999 0', '-SOLUTION/APRIORI',
                '+SOLUTION/ESTIMATE']
        for index, (axis, value) in enumerate(zip('XYZ', [6378137, 0, 0])):
            rows.append(f' {index+1} STA{axis} ABCD A 1 24:199:43200 m 2 {value} 0')
        if velocity:
            for index, axis in enumerate('XYZ'):
                rows.append(f' {index+4} VEL{axis} ABCD A 1 24:199:43200 m/y 2 1 0')
        rows.append('-SOLUTION/ESTIMATE')
        return '\n'.join(rows)

    def test_estimate_not_apriori_and_no_fake_velocity(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'fixture.snx'; path.write_text(self.fixture())
            groups, intervals = audit.sinex_positions(path)
            reference = audit.reference_for('ABCD', datetime(2024, 7, 17, 2), groups, intervals)
            self.assertEqual(reference['xyz_m'][0], 6378137)
            self.assertFalse(reference['velocity_propagated'])
            self.assertEqual(reference['epoch_difference_seconds'], -36000)
            self.assertIn('UNPROPAGATED', reference['epoch_scope'])

    def test_velocity_propagation(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'fixture.snx'; path.write_text(self.fixture(True))
            groups, intervals = audit.sinex_positions(path)
            reference = audit.reference_for('ABCD', datetime(2024, 7, 18, 12), groups, intervals)
            self.assertTrue(reference['velocity_propagated'])
            self.assertAlmostEqual(reference['xyz_m'][1], 1 / 365.25)

    def test_ambiguous_solution_not_silently_selected(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'fixture.snx'
            path.write_text(self.fixture() + '\n' + self.fixture().replace(' A 1 ', ' A 2 '))
            groups, intervals = audit.sinex_positions(path)
            reference = audit.reference_for('ABCD', datetime(2024, 7, 17, 2), groups, intervals)
            self.assertEqual(reference['status'], 'AMBIGUOUS_REFERENCE_SOLUTION')

    def test_missing_station_preserved_and_prior_not_independent(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'fixture.snx'; path.write_text(self.fixture())
            inventory = {'producer_sha256': 'producer', 'bundle_sha256': 'bundle',
                         'experiment_mode': audit.EXPECTED_MODE,
                         'checkpoint_epoch': '2024-07-17 02:00:00', 'coordinates': {},
                         'full_x_finite': True, 'full_P_finite': True,
                         'receiver_position_rate_count': 0, 'state_counts_by_system_and_type': {}}
            plan = {'binary_sha256': 'producer', 'pinned_input_contract':
                    {'observations': [{'station': 'ABCD'}, {'station': 'MISS'}]}}
            report = audit.audit(inventory, plan, path)
            self.assertEqual(report['expected_station_count'], 2)
            self.assertFalse(report['all_expected_XYZ_positive_finite_coordinate_covariance'])
            self.assertEqual(report['station_results']['MISS']['status'], 'NO_REC_POS_STATE')
            self.assertIn('NOT_INDEPENDENT', report['reference_claim'])
            json.dumps(report, allow_nan=False)

    def test_producer_mismatch_rejected(self):
        with self.assertRaisesRegex(ValueError, 'hashes differ'):
            audit.audit({'producer_sha256': 'a'}, {'binary_sha256': 'b'}, Path('unused'))

    def test_wrong_checkpoint_mode_rejected(self):
        with self.assertRaisesRegex(ValueError, 'target G/E/C FLOAT'):
            audit.audit({'producer_sha256': 'a', 'experiment_mode': 'E29_GPS_L1C_L2W_ZHANG_FULL_RANK'},
                        {'binary_sha256': 'a'}, Path('unused'))


if __name__ == '__main__':
    unittest.main()
