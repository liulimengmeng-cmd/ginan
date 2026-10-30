#!/usr/bin/env python3
"""Sample the isolated PEA run; stop it after sustained memory exhaustion."""

import datetime as dt
import json
import os
from pathlib import Path
import signal
import sys
import time


BINARY = b"/home/rx/GINAN/build-r51-p0p4-20260927/bin/pea\0"
RUN_ROOT = Path("/home/rx/GINAN/local_timing/r51_20260928_epoch_only_ar1h_01")


def kib_field(path: Path, name: str) -> int:
    for line in path.read_text().splitlines():
        if line.startswith(name + ":"):
            return int(line.split()[1])
    raise RuntimeError(f"Missing {name} in {path}")


def alive(pid: int) -> bool:
    try:
        return (Path("/proc") / str(pid) / "cmdline").read_bytes().startswith(BINARY)
    except FileNotFoundError:
        return False


def main() -> int:
    pid = int(sys.argv[1])
    if not RUN_ROOT.is_dir() or not alive(pid):
        raise RuntimeError("Expected isolated PEA process is not running")
    proc_status = Path("/proc") / str(pid) / "status"
    low_memory = 0
    high_swap = 0
    with (RUN_ROOT / "process_memory.jsonl").open("x", buffering=1) as output:
        while alive(pid):
            swap = kib_field(proc_status, "VmSwap")
            hwm = kib_field(proc_status, "VmHWM")
            rss = kib_field(proc_status, "VmRSS")
            available = kib_field(Path("/proc/meminfo"), "MemAvailable")
            low_memory = low_memory + 1 if available <= 1024 * 1024 else 0
            high_swap = high_swap + 1 if swap >= 2 * 1024 * 1024 else 0
            record = {
                "time": dt.datetime.now(dt.timezone.utc).isoformat(),
                "pid": pid,
                "rss_kib": rss,
                "hwm_kib": hwm,
                "swap_kib": swap,
                "mem_available_kib": available,
            }
            if low_memory >= 12 or high_swap >= 12:
                record["stop_reason"] = (
                    "SUSTAINED_AVAILABLE_MEMORY_UNDER_1_GIB"
                    if low_memory >= 12 else "SUSTAINED_SWAP_OVER_2_GIB"
                )
                output.write(json.dumps(record) + "\n")
                if alive(pid):
                    os.kill(pid, signal.SIGTERM)
                return 2
            output.write(json.dumps(record) + "\n")
            time.sleep(10)
        output.write(json.dumps({"ended": "target process exited"}) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
