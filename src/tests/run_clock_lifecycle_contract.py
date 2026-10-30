"""Actual PEA-linked controlled clock lifecycle; never launches a PPP job."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("build", type=Path)
parser.add_argument("--output", type=Path)
parser.add_argument("--syntax-only", action="store_true")
parser.add_argument("--build-only", action="store_true")
parser.add_argument("--observe-only", action="store_true")
args = parser.parse_args()
root = Path(__file__).resolve().parents[2]
build = args.build.resolve()
cache = (build / "CMakeCache.txt").read_text()
assert "CMAKE_HOME_DIRECTORY:INTERNAL=" + str(root / "src") in cache, "wrong checkout's PEA objects"
output = args.output.resolve() if args.output else build / "clock_lifecycle_contract"
output.mkdir(parents=True, exist_ok=True)
flags = {}
for line in (build / "cpp/CMakeFiles/pea.dir/flags.make").read_text().splitlines():
    if " = " in line:
        key, value = line.split(" = ", 1)
        flags[key] = shlex.split(value)
link = shlex.split((build / "cpp/CMakeFiles/pea.dir/link.txt").read_text())
prefix = [link[0], *flags["CXX_DEFINES"], *flags["CXX_INCLUDES"], *flags["CXX_FLAGS"], "-O0", "-std=c++20"]
driver = root / "src/tests/unit/multi_gnss_clock_lifecycle.cpp"
driver_hash = hashlib.sha256(driver.read_bytes()).hexdigest()
if args.syntax_only:
    raise SystemExit(subprocess.call([*prefix, "-fsyntax-only", str(driver)]))
objects = [build / "cpp" / item for item in link if item.endswith(".o")]
def object_stamps():
    return {str(path.relative_to(build)): [path.stat().st_size, path.stat().st_mtime_ns] for path in objects}
before_objects = object_stamps()
boundary_headers = [root / "src/cpp/common/receiver.hpp", root / "src/cpp/common/algebra.hpp",
                    root / "src/cpp/common/zhangSystemClockDatum.hpp",
                    root / "src/cpp/common/zhangEffectiveClockDatum.hpp",
                    root / "src/cpp/common/zhangClockLifecycle.hpp"]
before_headers = {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
                  for path in boundary_headers}
with (output / "build.log").open("w") as log:
    subprocess.run([*prefix, "-c", str(driver), "-o", str(output / "driver.o")],
                   stdout=log, stderr=subprocess.STDOUT, check=True)
    subprocess.run(["objcopy", "--redefine-sym", "main=pea_original_main",
        str(build / "cpp/CMakeFiles/pea.dir/pea/main.cpp.o"), str(output / "original_main.o")],
        stdout=log, stderr=subprocess.STDOUT, check=True)
    link[link.index("CMakeFiles/pea.dir/pea/main.cpp.o")] = str(output / "original_main.o")
    link.insert(1, str(output / "driver.o"))
    executable = output / "clock_lifecycle_contract"
    link[link.index("-o") + 1] = str(executable)
    subprocess.run(link, cwd=build / "cpp", stdout=log, stderr=subprocess.STDOUT, check=True)
assert before_objects == object_stamps(), "PEA objects changed during contract compilation/link; retry after the production build"
assert driver_hash == hashlib.sha256(driver.read_bytes()).hexdigest(), "contract driver changed during compilation; rerun the stable source"
assert before_headers == {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
                          for path in boundary_headers}, "clock ABI headers changed while compiling the contract"
(output / "linked_objects.json").write_text(json.dumps(before_objects, indent=2) + "\n")
if args.build_only:
    print(executable)
    raise SystemExit(0)
command = [str(executable), str(output / "runtime.TRACE")]
if args.observe_only:
    command.append("--observe-only")
with (output / "runtime.log").open("w") as log:
    result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
        env={**os.environ, "ZHANG_R51_RATIO_ONLY": "1", "ZHANG_R51_RATIO_THRESHOLD": "2",
             "ZHANG_R51_EPOCH_ONLY_AR": "1", "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1",
             "ZHANG_CLOCK_LIFECYCLE_DIAGNOSTICS": "1"})
report = {"exit_code": result.returncode, "observe_only": args.observe_only,
          "observe_only_is_acceptance": False, "real_ppp_run": False,
          "driver_sha256": driver_hash,
          "production_binary_sha256": hashlib.sha256((root / "bin/pea").read_bytes()).hexdigest(),
          "boundary_header_sha256": before_headers,
          "linked_objects_unchanged_during_build": True,
          "contract_executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
          "production_objects": str(build), "trace": str(output / "runtime.TRACE")}
(output / "result.json").write_text(json.dumps(report, indent=2) + "\n")
for line in (output / "runtime.log").read_text().splitlines():
    if line.startswith(("CLOCK_LIFECYCLE ", "PASS ", "FAIL ", "OBSERVED_FAILURE ")):
        print(line)
print(json.dumps(report))
raise SystemExit(result.returncode)
