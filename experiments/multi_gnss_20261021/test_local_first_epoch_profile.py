"""Guard-only fixtures; no GNSS numerical process is launched by these tests."""
import json
import os
from pathlib import Path
import signal
import sys
import tempfile
import unittest
from unittest.mock import patch

import run_local_first_epoch_profile as profile


class LocalProfileContracts(unittest.TestCase):
    def test_profile_flags_cannot_override_physics_or_duplicate(self):
        self.assertEqual(profile.profile_environment(["ZHANG_FLOAT_PROFILE=1"]), {"ZHANG_FLOAT_PROFILE": "1"})
        for value in (["ZHANG_R51_RATIO_THRESHOLD=1"], ["OMP_NUM_THREADS=100"],
                      ["ZHANG_FLOAT_PROFILE=1", "ZHANG_FLOAT_PROFILE=0"], ["ZHANG_FLOAT_PROFILE=\n"]):
            with self.assertRaises(ValueError):
                profile.profile_environment(value)

    def test_resource_boundaries_are_local_not_server_limits(self):
        sample = {"owned_tree_rss_bytes": 12 * profile.GIB - 1,
                  "host_available_memory_bytes": 4 * profile.GIB,
                  "output_free_disk_bytes": 20 * profile.GIB}
        self.assertEqual(profile.violations(sample), [])
        sample.update(owned_tree_rss_bytes=12 * profile.GIB,
            host_available_memory_bytes=4 * profile.GIB - 1, output_free_disk_bytes=20 * profile.GIB - 1)
        self.assertEqual(profile.violations(sample), ["OWNED_RSS_LIMIT", "HOST_MEMAVAILABLE_RESERVE", "OUTPUT_DISK_RESERVE"])

    def test_pid_reuse_or_wrong_session_is_never_signalled(self):
        saved = {"pid": 123, "start_ticks": 1000, "session": 123}
        for changed in (dict(saved, start_ticks=1001), dict(saved, session=999), None):
            with patch.object(profile, "proc_record", return_value=changed), \
                 patch.object(profile.os, "kill") as kill, patch.object(profile.time, "sleep"):
                self.assertEqual(profile.stop_owned([saved]), [])
                kill.assert_not_called()

    def test_identity_verified_term_then_kill_budget(self):
        saved = {"pid": 123, "start_ticks": 1000, "session": 123}
        with patch.object(profile, "proc_record", return_value=saved), \
             patch.object(profile.os, "kill") as kill, patch.object(profile.time, "sleep") as sleep:
            events = profile.stop_owned([saved])
        self.assertEqual([call.args[1] for call in kill.call_args_list], [signal.SIGTERM, signal.SIGKILL])
        self.assertEqual([call.args[0] for call in sleep.call_args_list], [1, 1])
        self.assertEqual([event["signal"] for event in events], ["SIGTERM", "SIGKILL"])

    @unittest.skipUnless(sys.platform == "linux", "Native Linux pidfd contract")
    def test_rss_capture_signals_only_verified_owned_gdb_and_waits_one_second(self):
        gdb = {"pid": 123, "start_ticks": 1000, "session": 123, "command": "gdb"}
        pea = dict(gdb, pid=124, command="pea")
        sample = {"host_available_memory_bytes": 4 * profile.GIB}
        with patch.object(profile, "proc_record", return_value=gdb), \
             patch.object(profile.os, "pidfd_open", return_value=7) as opened, \
             patch.object(profile.os, "close") as close, \
             patch.object(profile.signal, "pidfd_send_signal") as send, \
             patch.object(profile.time, "sleep") as sleep:
            events = profile.capture_rss_stop_stack([pea, gdb], sample, ["OWNED_RSS_LIMIT"])
        opened.assert_called_once_with(123)
        send.assert_called_once_with(7, signal.SIGINT, None, 0)
        close.assert_called_once_with(7)
        sleep.assert_called_once_with(1)
        self.assertEqual(events[0]["maximum_wait_seconds"], 1)
        self.assertEqual(events[0]["reason"], "BOUNDED_RSS_STOP_STACK_CAPTURE")

    @unittest.skipUnless(sys.platform == "linux", "Native Linux pidfd contract")
    def test_no_stack_delay_for_low_memory_non_rss_or_changed_identity(self):
        saved = {"pid": 123, "start_ticks": 1000, "session": 123, "command": "gdb"}
        with patch.object(profile.os, "pidfd_open", return_value=7) as opened, \
             patch.object(profile.os, "close"), \
             patch.object(profile.signal, "pidfd_send_signal") as send, \
             patch.object(profile.time, "sleep") as sleep, \
             patch.object(profile, "proc_record", return_value=dict(saved, start_ticks=1001)):
            for available, reasons in ((4 * profile.GIB - 1, ["OWNED_RSS_LIMIT"]),
                    (4 * profile.GIB, ["OUTPUT_DISK_RESERVE"]),
                    (4 * profile.GIB, ["OWNED_RSS_LIMIT"])):
                self.assertEqual(profile.capture_rss_stop_stack([saved],
                    {"host_available_memory_bytes": available}, reasons), [])
        opened.assert_called_once_with(123)
        send.assert_not_called()
        sleep.assert_not_called()

    @unittest.skipIf(profile.resource is None, "Linux resource module")
    def test_hard_as_and_core_limit_not_a_soft_rss_claim(self):
        with patch.object(profile.resource, "setrlimit") as limit:
            profile.hard_address_space_limit()
        self.assertEqual(limit.call_args_list[0].args,
            (profile.resource.RLIMIT_AS, (16 * profile.GIB, 16 * profile.GIB)))
        self.assertEqual(limit.call_args_list[1].args, (profile.resource.RLIMIT_CORE, (0, 0)))

    def test_gdb_is_parent_of_new_explicit_solver_not_attach(self):
        with patch.object(profile.network, "command", return_value=["/fake/pea", "--fixture"]):
            cmd = profile.profiler_command({"binary": "/fake/pea", "configs": [], "pinned_input_contract": {}}, True)
        self.assertIn("/usr/bin/gdb", cmd)
        self.assertIn("thread apply all bt 16", cmd)
        self.assertNotIn("attach", cmd)
        self.assertEqual(cmd[-3:], ["--args", "/fake/pea", "--fixture"])

    def test_duplicate_attempt_guard_precedes_all_execution(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            plan = {"state": "PREPARED_NOT_RUN_LOCAL_APPROVAL_REQUIRED", "run_root": str(root)}
            (root / "prepared_profile_plan.json").write_text(json.dumps(plan))
            (root / "stage_plan.json").write_text(json.dumps({"state": "PROFILE_RESOURCE_STOP"}))
            with patch.object(profile.subprocess, "Popen") as popen:
                with self.assertRaisesRegex(RuntimeError, "already attempted"):
                    profile.execute(plan, {}, root / "not_read_approval.json")
            popen.assert_not_called()

    @unittest.skipUnless(sys.platform == "linux", "Native Linux owned dummy process test")
    def test_monitor_discovery_failure_cleans_owned_dummy_and_records_failure(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            plan = {"run_root": str(root), "working_directory": str(root),
                "command": [sys.executable, "-c", "import time;time.sleep(30)"],
                "timeout_seconds": 2, "svd_profile_deadline_seconds": 15}
            with patch.object(profile, "owned_tree", return_value=[]):
                result = profile.supervise(plan, dict(os.environ))
            self.assertEqual(result["state"], "PROFILE_GUARD_FAILED")
            self.assertTrue(result["child_reaped"])
            self.assertEqual(result["owned_processes_still_live"], [])
            self.assertIn("monitoring unavailable", result["guard_failure_reason"])

    @unittest.skipUnless(sys.platform == "linux", "Native Linux owned dummy process test")
    def test_runtime_io_failure_cleans_owned_dummy_and_records_failure(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            plan = {"run_root": str(root), "working_directory": str(root),
                "command": [sys.executable, "-c", "import time;time.sleep(30)"],
                "timeout_seconds": 2, "svd_profile_deadline_seconds": 15}
            with patch.object(profile, "linux_memory", side_effect=OSError("INJECTED_MEMINFO_IO_FAILURE")):
                result = profile.supervise(plan, dict(os.environ))
            self.assertEqual(result["state"], "PROFILE_GUARD_FAILED")
            self.assertTrue(result["child_reaped"])
            self.assertEqual(result["owned_processes_still_live"], [])
            self.assertIn("INJECTED_MEMINFO_IO_FAILURE", result["guard_failure_reason"])

    @unittest.skipUnless(sys.platform == "linux", "Native Linux owned dummy session test")
    def test_initial_identity_unavailable_cleans_new_session_and_its_dummy_child(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            # Like a GDB inferior, the descendant has a different PGID but
            # remains in the newly owned SID. A killpg(root) would miss it.
            code = "import os,subprocess,sys,time; p=subprocess.Popen([sys.executable,'-c','import time;time.sleep(30)'],preexec_fn=os.setpgrp); print(p.pid,flush=True); time.sleep(30)"
            plan = {"run_root": str(root), "working_directory": str(root),
                "command": [sys.executable, "-c", code], "timeout_seconds": 2,
                "svd_profile_deadline_seconds": 15}
            real_record = profile.proc_record
            def unavailable(pid):
                # Allow the dummy to fork; this is test injection only.
                import time
                time.sleep(.15)
                return None
            with patch.object(profile, "proc_record", side_effect=unavailable):
                result = profile.supervise(plan, dict(os.environ))
            self.assertEqual(result["state"], "PROFILE_GUARD_FAILED")
            self.assertTrue(result["child_reaped"])
            self.assertTrue(any(row.get("ownership_basis", "").startswith("DIRECT_UNREAPED")
                                for row in result["identity_verified_signal_events"]))
            descendant = int((root / "runner.log").read_text().strip())
            self.assertTrue(any(row.get("pid") == descendant and row.get("signal") == "SIGTERM"
                                for row in result["identity_verified_signal_events"]))
            actual = real_record(descendant)
            self.assertTrue(actual is None or actual["state"] == "Z")


if __name__ == "__main__":
    unittest.main()
