"""Isolated one-committed-epoch probe; expected intentional exit code is 130.

Run as WSL root: python3 SCRIPT BINARY EVIDENCE_DIRECTORY UNIQUE_SUFFIX
Preserves the authenticated checkpoint/configuration. The existing launcher
binds only new output directories inside a private mount namespace. One
SIGINT after the first COMMIT asks PEA to exit at an epoch boundary.
"""
from pathlib import Path
import atexit
import hashlib
import json
import os
import re
import shutil
import signal
import pwd
import subprocess
import sys
import time


def sha256(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()

assert os.geteuid() == 0, "WSL root required for private bind mounts"
source, evidence = map(lambda s: Path(s).resolve(), sys.argv[1:3])
suffix = sys.argv[3]
assert re.fullmatch(r"[a-zA-Z0-9_-]+", suffix)
assert source.is_file() and source.name == "pea"
assert evidence.parent in (Path("/home/rx/GINAN/ar_optimization_20261019"),
                          Path("/home/rx/GINAN/product_closure_20261019"),
                          Path("/home/rx/GINAN/wl_l1_float_20261020"))
assert not evidence.exists(), "Refuse to replace an existing probe"
run = Path("/home/rx/GINAN/local_timing") / (
    "r51_20261001_osb_branch_repair_ar1h_" + suffix)
assert not run.exists()
assert subprocess.run(["pgrep", "-x", "pea"], stdout=subprocess.DEVNULL).returncode != 0, \
    "Another solver is running"
evidence.mkdir()
if evidence.parent.name=="wl_l1_float_20261020":
    condition_inputs=evidence/"condition_inputs"
    condition_inputs.mkdir()
    solver_account=pwd.getpwnam("rx")
    os.chown(condition_inputs,solver_account.pw_uid,solver_account.pw_gid)
    subprocess.run(["runuser","-u","rx","--","test","-w",str(condition_inputs)],check=True)
    os.environ["ZHANG_R51_CONDITION_REPLAY_DIRECTORY"]=str(condition_inputs)
    os.environ["ZHANG_R51_BLOCK_CONDITIONER"]="incremental"
    os.environ["ZHANG_R51_FLOAT_MATRIX_MODE"]="symmetric"
    exact_inputs=evidence/"exact_inputs"
    exact_inputs.mkdir()
    os.chown(exact_inputs,solver_account.pw_uid,solver_account.pw_gid)
    subprocess.run(["runuser","-u","rx","--","test","-w",str(exact_inputs)],check=True)
    os.environ["ZHANG_R51_EXACT_REPLAY_DIRECTORY"]=str(exact_inputs)
if evidence.parent.name=="product_closure_20261019":
    # Opt-in exact inputs can be replayed against a frozen old backend without
    # paying RINEX loading/checkpoint hashing for every microbenchmark.
    exact_inputs=evidence/"exact_inputs"
    exact_inputs.mkdir()
    # The namespace is created by root, but the authenticated launcher drops
    # the solver to rx. Only this NEW diagnostic output directory is writable
    # by that account; do not chown the binary, checkpoint or input tree.
    solver_account=pwd.getpwnam("rx")
    os.chown(exact_inputs,solver_account.pw_uid,solver_account.pw_gid)
    subprocess.run(["runuser","-u","rx","--","test","-w",str(exact_inputs)],check=True)
    os.environ["ZHANG_R51_EXACT_REPLAY_DIRECTORY"]=str(exact_inputs)
binary = evidence / "pea"
shutil.copy2(source, binary)
expected_hash = sha256(binary)
launcher = Path("/home/rx/GINAN/frozen-osb-ratio2-domain-budget-fix-20261001-v5/launcher.sh")
qc_config = Path("/home/rx/GINAN/frozen-cycle-slip-20261018-qualified-gf-pinned-v1/config")
config_target = Path("/home/rx/GINAN/frozen-r51-satellite-fix-20260918/config")
assert launcher.is_file() and qc_config.is_dir() and config_target.is_dir()
checkpoint = Path("/mnt/d/GINAN_R20/zhang_r51_ratio3_float2h_ar1h_2024199_180_20260917_checkpoints_source_before_ratioonly/E29-20240717015930-e240-rb06f1454a2e1/checkpoint.bundle")
checkpoint_hash = sha256(checkpoint)
assert checkpoint_hash == "ab397e655774da5c0a9632becd08fe171508b93ffeaf2dd7b04de899bae2182a", \
    "Original FLOAT checkpoint bytes changed"
events = []
manifest = dict(binary=str(binary), binary_sha256=expected_hash,
                run_root=str(run), maximum_seconds=1500, expected_returncode=128+signal.SIGINT,
                scope="one committed real epoch; intentional interruption, not a completed hour",
                checkpoint=str(checkpoint), checkpoint_sha256=checkpoint_hash,
                launcher_sha256=sha256(launcher),
                config_sha256={str(p.relative_to(qc_config)): sha256(p)
                               for p in sorted(qc_config.rglob("*.yaml"))},
                ratio=2, joint_nis=False, epoch_only=True, events=events)


def save():
    (evidence / "probe.json").write_text(json.dumps(manifest, indent=2))


def descendants(pid):
    result = []
    pending = [pid]
    for current in pending:
        result.append(current)
        try:
            children = Path(f"/proc/{current}/task/{current}/children").read_text()
            pending.extend(map(int, children.split()))
        except FileNotFoundError:
            pass
    return result


def owned_solver(pid):
    try:
        return Path(f"/proc/{pid}/exe").resolve() == binary
    except OSError:
        return False


save()
started = time.monotonic()
with (evidence / "launch.log").open("w") as log:
    process = subprocess.Popen([
        "timeout", "--signal=KILL", "26m", "unshare", "--mount", "--propagation", "private", "--", "bash", "-c",
        'set -euo pipefail; ulimit -c 0; mount --bind "$1" "$2"; '
        'export ZHANG_R51_EVIDENCE_SLIP_QC=1 ZHANG_R51_SLIP_RAW_WARMUP=1; '
        'exec bash "$3" "$4" "$5"',
        "probe", str(qc_config), str(config_target), str(launcher), str(binary), str(run)],
        stdout=log, stderr=subprocess.STDOUT)
    solver = None

    def cleanup_on_exception():
        # Only this probe's hash-verified executable, never all PEA processes.
        if solver is not None and owned_solver(solver):
            os.kill(solver, signal.SIGKILL)
        if process.poll() is None:
            process.terminate()
            process.wait(timeout=45)

    atexit.register(cleanup_on_exception)
    interrupted = False
    timed_out = False
    while process.poll() is None:
        if solver is None:
            candidates = [pid for pid in descendants(process.pid) if owned_solver(pid)]
            assert len(candidates) <= 1, "Multiple probe solvers"
            if candidates:
                solver = candidates[0]
                assert sha256(Path(f"/proc/{solver}/exe")) == expected_hash
                events.append(dict(event="verified_solver", pid=solver,
                                   elapsed=time.monotonic() - started))
                save()
                print(json.dumps(events[-1]), flush=True)
        commits = sorted((run / "outputs/zhang_internal_products.csv.epochs").glob("*/COMMIT.json"))
        ar_commits = []
        for path in commits:
            try:
                if json.loads(path.read_text()).get("ar") is True:
                    ar_commits.append(path)
            except (FileNotFoundError, json.JSONDecodeError):
                continue  # Wait for a complete atomic publication marker.
        if ar_commits and not interrupted and solver is not None and owned_solver(solver):
            os.kill(solver, signal.SIGINT)
            interrupted = True
            events.append(dict(event="single_SIGINT_after_AR_COMMIT", commit=str(ar_commits[0]),
                               elapsed=time.monotonic() - started))
            save()
            print(json.dumps(events[-1]), flush=True)
        if time.monotonic() - started > 1500:
            timed_out = True
            if solver is not None and owned_solver(solver):
                os.kill(solver, signal.SIGKILL)
            else:
                process.terminate()
            events.append(dict(event="probe_deadline", elapsed=time.monotonic() - started))
            break
        time.sleep(.2)
    returncode = process.wait(timeout=45)

commits = sorted((run / "outputs/zhang_internal_products.csv.epochs").glob("*/COMMIT.json"))
commit_records = [json.loads(path.read_text()) for path in commits]
manifest.update(returncode=returncode, timed_out=timed_out,
                elapsed_seconds=time.monotonic() - started,
                committed_epochs=[str(path) for path in commits],
                committed_epoch_times=sorted({row["epoch"] for row in commit_records}),
                ar_committed_bundle_count=sum(row.get("ar") is True for row in commit_records),
                checkpoint_bundles=[str(p) for p in sorted((run/"checkpoints").glob("*/checkpoint.bundle"))
                                    if (p.parent/"checkpoint_manifest.json").is_file()],
                solver_still_running=solver is not None and owned_solver(solver))
save()
print(json.dumps(manifest, indent=2), flush=True)
assert not timed_out and not manifest["solver_still_running"], "Probe did not exit cleanly"
assert returncode == manifest["expected_returncode"] and manifest["ar_committed_bundle_count"] > 0, \
    "Expected one AR epoch followed by intentional exit 130"
assert len(manifest["committed_epoch_times"]) == 1 and manifest["checkpoint_bundles"], \
    "Probe must finish one epoch and preserve its authenticated checkpoint"
