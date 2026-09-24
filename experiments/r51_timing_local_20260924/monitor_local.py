#!/usr/bin/env python3
"""Record one local PEA run and stop sustained memory exhaustion."""

from __future__ import annotations

import datetime as dt
import json
import os
from pathlib import Path
import signal
import sys
import time


EXPECTED_BINARY = b"/home/rx/GINAN/build-r51-timing-20260924/bin/pea"
SAMPLE_SECONDS = 10
CONSECUTIVE_LIMIT = 12
SWAP_LIMIT_KIB = 2 * 1024 * 1024
AVAILABLE_LIMIT_KIB = 1024 * 1024


def kib_field(path: Path, name: str) -> int:
    for line in path.read_text().splitlines():
        if line.startswith(name + ":"):
            return int(line.split()[1])
    raise RuntimeError(f"Missing {name} in {path}")


def target_alive(pid: int) -> bool:
    command = (Path("/proc") / str(pid) / "cmdline").read_bytes()
    return command.startswith(EXPECTED_BINARY + b"\0")


def main() -> int:
    pid = int(sys.argv[1])
    run_root = Path(sys.argv[2]).resolve()
    allowed_root = Path("/home/rx/GINAN/local_timing").resolve()
    if run_root.parent != allowed_root or not run_root.is_dir():
        raise RuntimeError(f"Unexpected run directory: {run_root}")
    output = run_root / "process_memory.jsonl"
    proc_status = Path("/proc") / str(pid) / "status"
    low_memory_samples = 0
    high_swap_samples = 0
    with output.open("x", buffering=1) as stream:
        while True:
            try:
                if not target_alive(pid):
                    stream.write(json.dumps({"ended": "target process exited"}) + "\n")
                    return 0
                swap_kib = kib_field(proc_status, "VmSwap")
                hwm_kib = kib_field(proc_status, "VmHWM")
                rss_kib = kib_field(proc_status, "VmRSS")
                available_kib = kib_field(Path("/proc/meminfo"), "MemAvailable")
            except (FileNotFoundError, ProcessLookupError):
                stream.write(json.dumps({"ended": "target process exited"}) + "\n")
                return 0
            high_swap_samples = high_swap_samples + 1 if swap_kib >= SWAP_LIMIT_KIB else 0
            low_memory_samples = low_memory_samples + 1 if available_kib <= AVAILABLE_LIMIT_KIB else 0
            record = {
                "time": dt.datetime.now(dt.timezone.utc).isoformat(),
                "pid": pid,
                "rss_kib": rss_kib,
                "hwm_kib": hwm_kib,
                "swap_kib": swap_kib,
                "mem_available_kib": available_kib,
            }
            if max(high_swap_samples, low_memory_samples) >= CONSECUTIVE_LIMIT:
                record["stop_reason"] = (
                    "SUSTAINED_SWAP_OVER_2_GIB" if high_swap_samples >= CONSECUTIVE_LIMIT
                    else "SUSTAINED_AVAILABLE_MEMORY_UNDER_1_GIB"
                )
                stream.write(json.dumps(record) + "\n")
                if target_alive(pid):
                    os.kill(pid, signal.SIGTERM)
                return 2
            stream.write(json.dumps(record) + "\n")
            time.sleep(SAMPLE_SECONDS)


if __name__ == "__main__":
    raise SystemExit(main())
