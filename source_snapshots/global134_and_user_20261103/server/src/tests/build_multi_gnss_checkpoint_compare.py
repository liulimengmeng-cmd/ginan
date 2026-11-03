"""Compile only the independent test reader; link existing PEA objects.

No production source is recompiled and no solver is launched.
Usage: BUILD_DIR OUTPUT_DIR [--compile-only] [--blas-library /absolute/libopenblas.so]
"""
from pathlib import Path
import argparse
import shlex
import subprocess
import sys

repo = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("build", type=Path)
parser.add_argument("output", type=Path)
parser.add_argument("--compile-only", action="store_true")
parser.add_argument("--blas-library", type=Path)
parser.add_argument("--source", type=Path, default=repo / "src/tests/unit/multi_gnss_checkpoint_compare.cpp")
parser.add_argument("--executable-name", default="multi_gnss_checkpoint_compare")
parser.add_argument("--linear-combo-object", type=Path)
args = parser.parse_args()
build, output = args.build.resolve(), args.output.resolve()
output.mkdir(parents=True, exist_ok=True)
flags = {}
for line in (build / "cpp/CMakeFiles/pea.dir/flags.make").read_text().splitlines():
    if " = " in line:
        name, value = line.split(" = ", 1)
        flags[name] = shlex.split(value)
command = shlex.split((build / "cpp/CMakeFiles/pea.dir/link.txt").read_text())
if args.blas_library is not None:
    library = args.blas_library.resolve(strict=True)
    selected = [index for index, value in enumerate(command) if "openblas" in Path(value).name.lower()]
    if not selected:
        raise RuntimeError("No explicit OpenBLAS link argument found; refusing an unreviewed fallback")
    for index in selected:
        command[index] = str(library)
    command.append("-Wl,-rpath," + str(library.parent))
obj = output / "multi_gnss_checkpoint_compare.o"
subprocess.run([command[0], *flags["CXX_DEFINES"], *flags["CXX_INCLUDES"], *flags["CXX_FLAGS"],
    "-O2", "-std=c++20", "-c", str(args.source.resolve(strict=True)), "-o", str(obj)], check=True)
if args.compile_only:
    print(obj, flush=True)
    raise SystemExit(0)
renamed = output / "pea_main_renamed.o"
subprocess.run(["objcopy", "--redefine-sym", "main=pea_original_main",
    str(build / "cpp/CMakeFiles/pea.dir/pea/main.cpp.o"), str(renamed)], check=True)
command[command.index("CMakeFiles/pea.dir/pea/main.cpp.o")] = str(renamed)
if args.linear_combo_object is not None:
    command[command.index("CMakeFiles/pea.dir/common/linearCombo.cpp.o")] = str(args.linear_combo_object.resolve(strict=True))
command.insert(1, str(obj))
if Path(args.executable_name).name != args.executable_name:
    raise RuntimeError("Executable name must be a basename")
executable = output / args.executable_name
command[command.index("-o") + 1] = str(executable)
subprocess.run(command, cwd=build / "cpp", check=True)
print(executable, flush=True)
