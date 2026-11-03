"""Link/run the physical-consequence oracle against the actual PEA objects.

Usage: python3 src/tests/run_r51_physical_consequence_production.py <build-dir>
This does not launch the solver, alter existing objects, or replace pea.
The caller must finish the candidate PEA build first.
"""
from pathlib import Path
import os
import shlex
import subprocess
import sys


def main():
    root = Path(__file__).resolve().parents[2]
    build = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else root / "build"
    output = build / "physical_consequence_production_contract"
    output.mkdir(exist_ok=True)
    flags = {}
    for line in (build / "cpp/CMakeFiles/pea.dir/flags.make").read_text().splitlines():
        if " = " in line:
            key, value = line.split(" = ", 1)
            flags[key] = shlex.split(value)
    command = shlex.split((build / "cpp/CMakeFiles/pea.dir/link.txt").read_text())
    driver = output / "driver.o"
    subprocess.run([
        command[0], *flags["CXX_DEFINES"], *flags["CXX_INCLUDES"], *flags["CXX_FLAGS"],
        "-O1", "-fopenmp", "-std=c++20", "-c",
        str(root / "src/tests/unit/r51_physical_consequence_production.cpp"),
        "-o", str(driver)], check=True)
    original_main = output / "original_main.o"
    subprocess.run([
        "objcopy", "--redefine-sym", "main=pea_original_main",
        str(build / "cpp/CMakeFiles/pea.dir/pea/main.cpp.o"),
        str(original_main)], check=True)
    command[command.index("CMakeFiles/pea.dir/pea/main.cpp.o")] = str(original_main)
    command.insert(1, str(driver))
    executable = output / "physical_consequence_production_contract"
    command[command.index("-o") + 1] = str(executable)
    subprocess.run(command, cwd=build / "cpp", check=True)
    subprocess.run([str(executable)], check=True, env={
        **os.environ, "ZHANG_R51_ENABLE": "1", "ZHANG_R51_RATIO_ONLY": "1",
        "ZHANG_R51_RATIO_THRESHOLD": "2", "ZHANG_R51_EPOCH_ONLY_AR": "1",
        "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"})


if __name__ == "__main__":
    main()
