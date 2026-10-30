"""Server supervisor guard contracts; no GNSS numerical process."""
import json
from pathlib import Path
import signal
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch
import run_server_first_epoch_profile as server


class ServerProfileContracts(unittest.TestCase):
    def test_server_limits_are_explicit_and_not_local_limits(self):
        self.assertEqual(server.SERVER_LIMITS["maximum_address_space_bytes"], 64 * server.GIB)
        self.assertEqual(server.SERVER_LIMITS["maximum_owned_tree_rss_bytes"], 40 * server.GIB)
        self.assertEqual(server.SERVER_LIMITS["minimum_host_available_memory_bytes"], 16 * server.GIB)
        self.assertEqual(server.SERVER_LIMITS["sample_interval_seconds"], 1)

    def test_portable_loader_prefix_and_main_thread_stack_precede_all_threads(self):
        prefix = server.debugger_prefix(Path("/private/gdb"))
        plan = {"debugger_prefix": prefix, "binary": "/private/pea", "configs": [], "pinned_input_contract": {}}
        network = SimpleNamespace(command=lambda *args: ["/private/pea", "--fixture"])
        cmd = server.numerical_command(plan, network)
        self.assertEqual(cmd[2:2 + len(prefix)], prefix)
        self.assertLess(cmd.index("thread 1"), cmd.index("bt 32"))
        self.assertLess(cmd.index("bt 32"), cmd.index("thread apply all bt 8"))
        self.assertIn("-iex", cmd)
        self.assertIn("set debuginfod enabled off", cmd)
        self.assertNotIn("attach", cmd)

    def test_only_exact_argv_birth_and_session_can_receive_debugger_interrupt(self):
        saved = {"pid": 123, "start_ticks": 100, "session": 120, "command": "ld-linux-x86-64"}
        prefix = ["/loader", "--library-path", "/private/lib", "/private/bin/gdb"]
        for actual, argv in ((dict(saved, start_ticks=101), prefix),
                             (dict(saved, session=999), prefix), (saved, ["/loader", "/other/gdb"])):
            profile = SimpleNamespace(proc_record=lambda pid, actual=actual: actual)
            with patch.object(server, "read_argv", return_value=argv), \
                 patch.object(server.os, "pidfd_open", return_value=7), \
                 patch.object(server.os, "close"), \
                 patch.object(server.signal, "pidfd_send_signal") as sent:
                with self.assertRaises(RuntimeError):
                    server.interrupt_debugger([saved], prefix, profile)
                sent.assert_not_called()

    def test_loader_backed_gdb_exact_identity_receives_pidfd_sigint(self):
        saved = {"pid": 123, "start_ticks": 100, "session": 120, "command": "ld-linux-x86-64"}
        prefix = ["/loader", "--library-path", "/private/lib", "/private/bin/gdb"]
        profile = SimpleNamespace(proc_record=lambda pid: saved)
        with patch.object(server, "read_argv", return_value=prefix + ["--batch"]), \
             patch.object(server.os, "pidfd_open", return_value=7), \
             patch.object(server.os, "close"), \
             patch.object(server.signal, "pidfd_send_signal") as sent:
            event = server.interrupt_debugger([saved], prefix, profile)
        sent.assert_called_once_with(7, signal.SIGINT, None, 0)
        self.assertEqual(event["actual_command"], "ld-linux-x86-64")

    def test_completed_svd_never_requests_a_later_stage_svd_stack(self):
        end = 'ZHANG_FLOAT_PHASE_END phase=FLOAT_RANK_JACOBI_SVD epoch="2024-07-17 00:00:00"'
        self.assertEqual(server.sampling_action(end, 0, None, 16), "SVD_COMPLETED_BEFORE_STACK_SAMPLE")
        self.assertIsNone(server.sampling_action("", 0, None, 14.9))
        self.assertEqual(server.sampling_action("", 0, None, 15), "REQUEST_SVD_STACK")

    def test_semantic_overlay_diff_exposes_extra_threshold_change(self):
        old = {"processing_options": {"ratio": 2}, "outputs": {"outputs_root": "/old"}}
        new = {"processing_options": {"ratio": 1}, "outputs": {"outputs_root": "/new"}}
        changes = server.semantic_changes(old, new)
        self.assertEqual(changes[0]["path"], ["processing_options", "ratio"])
        self.assertEqual(changes[0]["new"], 1)

    def test_duplicate_attempt_rejected_before_approval_or_popen(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            plan = {"run_root": str(root)}
            (root / "prepared_server_profile_plan.json").write_text(json.dumps(plan))
            (root / "stage_plan.json").write_text(json.dumps({"state": "PROFILE_RESOURCE_STOP"}))
            with patch.object(server.subprocess, "Popen") as popen:
                with self.assertRaisesRegex(RuntimeError, "unused exact"):
                    server.execute(plan, root / "not_read.json", None, None, None)
                popen.assert_not_called()


if __name__ == "__main__":
    unittest.main()
