import unittest
import tempfile
from pathlib import Path
from experiment5_code_osb_sigma import perturb, C
from audit_experiment5_code_osb_sigma import screening


def line(sat='G01', signal='L1C', sigma=0):
    b = bytearray(b' '*103 + b'\r\n')
    for start, value in [(1,b'OSB '),(11,sat.encode()),(25,signal.encode()),
                         (65,b'ns '),(70,b'              12.3456'),
                         (92,f'{sigma:11.7f}'.encode())]:
        b[start:start+len(value)] = value
    return bytes(b)


class FixtureTests(unittest.TestCase):
    def setUp(self):
        self.target = b''.join(line(f'G{s:02}', sig) for s in range(1,33) for sig in ('L1C','L2W'))

    def test_only_phase_sigma_bytes_change(self):
        extra = b'* retain provenance\r\n' + line(signal='C1C',sigma=.006) + line(sat='E01')
        original = extra+self.target
        derived, changes = perturb(original,.003)
        self.assertEqual(derived[:len(extra)],extra)
        self.assertEqual(len(changes),64)
        for old,new in zip(original.splitlines(True),derived.splitlines(True)):
            self.assertEqual(old[:92],new[:92])
            self.assertEqual(old[103:],new[103:])
        self.assertAlmostEqual(changes[0]['effective_sigma_m'],.003,delta=C/1e9*0.5e-7)

    def test_rejects_nonzero_original_sigma(self):
        with self.assertRaises(AssertionError):
            perturb(self.target.replace(b'  0.0000000',b'  0.0010000',1),.01)

    def test_rejects_missing_records(self):
        with self.assertRaises(AssertionError):
            perturb(self.target.split(b'\r\n',1)[1],.01)

    def test_rejects_nonpositive_sigma(self):
        with self.assertRaises(AssertionError):
            perturb(self.target,0)

    def test_screening_is_scoped_to_ar_shadow_call(self):
        text = ('Postfit check failed measurement test\n'
                'Epoch 12 = test\n'
                'PPP_AR FEEDBACK_SHADOW_SUMMARY rows=1\n'
                'Postfit check failed measurement test\n'
                'PPP_AR FEEDBACK_SHADOW_REALISATION state_update_error_norm=0.5\n'
                'Postfit check failed measurement test\n'
                'Epoch 13 = test\n'
                'PPP_AR FEEDBACK_SHADOW_SUMMARY rows=1\n'
                'PPP_AR FEEDBACK_SHADOW_REALISATION state_update_error_norm=1e-9\n')
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'trace'
            p.write_text(text)
            result=screening(p)
        self.assertEqual(result['postfit_failure_events'],1)
        self.assertEqual(result['calls_with_screening_failure'],1)
        self.assertEqual(result['feedback_calls'],2)
        self.assertEqual(result['mixed_unit_norm_without_screening']['max'],1e-9)

if __name__ == '__main__':
    unittest.main()
