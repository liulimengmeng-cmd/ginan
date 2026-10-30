"""Candidate audit reader invariants; these tests do not run PEA."""
import json
from pathlib import Path
import tempfile
import unittest

import audit_candidate_replays as audit


class CandidateAuditTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_running_snapshot_never_promoted_to_completion(self):
        (self.root / "stage_plan.json").write_text(json.dumps({"state": "RUNNING"}))
        result = audit.audit(self.root)
        self.assertEqual(result["completion"], "NOT_ESTABLISHED")
        self.assertFalse(result["scientific_acceptance"])

    def test_joint_clock_formula_uses_explicit_system_value_not_selected_system(self):
        row = {"source_line": 77, "jsonl_line": 1, "event": {"phase": "REC_CLOCK_ADJUST_AFTER",
            "receiver": "NNOR", "epoch": "t1", "system": "GAL", "spp_system": "GPS",
            "spp_clock_m": 9999, "spp_system_clock_m": 64, "saved_offset_m": 9999,
            "saved_system_offset_m": 647, "predicted_clock_m": 712,
            "unrounded_adjustment_m": -1, "applied_adjustment_m": -1}}
        result = audit.offset_chain([row])["formula_checks"][0]
        self.assertEqual(result["system"], "GAL")
        self.assertEqual(result["independently_recomputed_adjustment_m"], -1)
        self.assertEqual(result["formula_residual_m"], 0)

    def test_two_point_raw_phase_change_is_not_accepted_continuity(self):
        for epoch, phase in (("t0", 2), ("t1", 65)):
            root = self.root / "outputs/zhang_internal_products.csv.epochs" / epoch
            root.mkdir(parents=True)
            (root / "COMMIT.json").write_text("{}")
            (root / "products.csv").write_text("solution,satellite,observable,phase_m\nFLOAT,E02,L1C," + str(phase) + "\n")
        result = audit.raw_phase_changes(self.root)
        self.assertEqual(result["maximum_absolute_difference_m"], 63)
        self.assertIn("pending convergence", result["interpretation"])
        self.assertNotIn("scientific_acceptance", result)


if __name__ == "__main__":
    unittest.main()
