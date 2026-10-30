"""Build only a checkpoint test DRIVER and link existing production objects.

No production C++ file is compiled and no solver is launched.
Usage: BUILD OUTPUT [real.bundle producer-binary-sha256 [--upgrade-real-schema5]]
"""
from pathlib import Path
import shlex
import subprocess
import sys

repo = Path(__file__).resolve().parents[2]
build = Path(sys.argv[1]).resolve()
output = Path(sys.argv[2]).resolve()
output.mkdir(parents=True, exist_ok=True)
flags = {}
for line in (build / "cpp/CMakeFiles/pea.dir/flags.make").read_text().splitlines():
    if " = " in line:
        key, value = line.split(" = ", 1)
        flags[key] = shlex.split(value)
command = shlex.split((build / "cpp/CMakeFiles/pea.dir/link.txt").read_text())
driver = output / "structural_checkpoint_driver.o"
subprocess.run([command[0], *flags["CXX_DEFINES"], *flags["CXX_INCLUDES"],
    *flags["CXX_FLAGS"], "-O2", "-std=c++20", "-c",
    str(repo / "src/tests/unit/r51_structural_checkpoint_roundtrip.cpp"),
    "-o", str(driver)], check=True)
renamed_main = output / "pea_main_renamed.o"
subprocess.run(["objcopy", "--redefine-sym", "main=pea_original_main",
    str(build / "cpp/CMakeFiles/pea.dir/pea/main.cpp.o"), str(renamed_main)], check=True)
command[command.index("CMakeFiles/pea.dir/pea/main.cpp.o")] = str(renamed_main)
command.insert(1, str(driver))
executable = output / "structural_checkpoint_roundtrip"
command[command.index("-o") + 1] = str(executable)
subprocess.run(command, cwd=build / "cpp", check=True)
print(executable, flush=True)
if len(sys.argv) > 3:
    assert len(sys.argv) in (5, 6), "Provide real bundle, producer SHA256, optional schema5 upgrade mode"
    subprocess.run([str(executable), *sys.argv[3:]], check=True)
