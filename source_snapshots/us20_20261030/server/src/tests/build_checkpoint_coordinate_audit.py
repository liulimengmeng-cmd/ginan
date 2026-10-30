"""Link a read-only checkpoint inventory driver to stable production objects."""
from pathlib import Path
import argparse
import shlex
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('build', type=Path)
parser.add_argument('output', type=Path)
parser.add_argument('--syntax-only', action='store_true')
parser.add_argument('--blas-library', type=Path)
parser.add_argument('--clock-metadata-only', action='store_true')
args = parser.parse_args()
repo = Path(__file__).resolve().parents[2]
build = args.build.resolve()
output = args.output.resolve()
output.mkdir(parents=True, exist_ok=True)
settings = {}
for line in (build / 'cpp/CMakeFiles/pea.dir/flags.make').read_text().splitlines():
    key, separator, value = line.partition(' = ')
    if separator:
        settings[key] = shlex.split(value)
link = shlex.split((build / 'cpp/CMakeFiles/pea.dir/link.txt').read_text())
if args.blas_library is not None:
    library = args.blas_library.resolve(strict=True)
    selected = [i for i, item in enumerate(link) if 'openblas' in Path(item).name.lower()]
    if not selected:
        raise RuntimeError('No explicit OpenBLAS link argument; refusing an unreviewed fallback')
    for index in selected:
        link[index] = str(library)
    link.append('-Wl,-rpath,' + str(library.parent))
driver_name = 'multi_gnss_checkpoint_clock_metadata_audit.cpp' if args.clock_metadata_only else 'multi_gnss_checkpoint_coordinate_audit.cpp'
executable_name = 'checkpoint_clock_metadata_inventory' if args.clock_metadata_only else 'checkpoint_coordinate_inventory'
driver = repo / 'src/tests/unit' / driver_name
prefix = [link[0], *settings['CXX_DEFINES'], *settings['CXX_INCLUDES'], *settings['CXX_FLAGS'], '-O0']
if args.syntax_only:
    raise SystemExit(subprocess.call([*prefix, '-fsyntax-only', str(driver)]))
objects = [build / 'cpp' / item for item in link if item.endswith('.o')]
def stamps():
    return [(str(path), path.stat().st_size, path.stat().st_mtime_ns) for path in objects]
before = stamps()
with (output / 'build.log').open('w') as log:
    obj = output / 'coordinate_inventory.o'
    subprocess.run([*prefix, '-c', str(driver), '-o', str(obj)], check=True,
                   stdout=log, stderr=subprocess.STDOUT)
    renamed = output / 'original_pea_main.o'
    subprocess.run(['objcopy', '--redefine-sym', 'main=pea_original_main',
                    str(build / 'cpp/CMakeFiles/pea.dir/pea/main.cpp.o'), str(renamed)],
                   check=True, stdout=log, stderr=subprocess.STDOUT)
    link[link.index('CMakeFiles/pea.dir/pea/main.cpp.o')] = str(renamed)
    link.insert(1, str(obj))
    link[link.index('-o') + 1] = str(output / executable_name)
    subprocess.run(link, cwd=build / 'cpp', stdout=log, stderr=subprocess.STDOUT, check=True)
assert before == stamps(), 'Production objects changed while linking inventory reader'
print(output / executable_name)
