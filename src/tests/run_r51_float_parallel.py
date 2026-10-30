"""Build a production-linked FLOAT benchmark; optional 'baseline' uses old objects.

Usage: python3 src/tests/run_r51_float_parallel.py BUILD OUTPUT [baseline|frozen] [ALGEBRA_OBJECT]
       [--blas-library /absolute/path/to/libopenblas.so]
The benchmark executable is emitted, not automatically run. PEA is not modified.
"""
import argparse
from pathlib import Path
import shlex
import subprocess


def profile_definitions(mode):
    # Retain the original positional baseline/frozen object-linking contract.
    result = [] if mode == 'baseline' else ['-DFLOAT_THREAD_PROFILE=1']
    if mode == 'frozen':
        result.append('-DFLOAT_FROZEN_PROFILE=1')
    return result


def compatible_blas_link(command, library):
    """Replace every explicit OpenBLAS archive/library, never unrelated flags.

    No implicit installed-library fallback: the caller selects and freezes the
    private generic library. Its parent is added to the executable's RPATH.
    """
    library = Path(library).resolve(strict=True)
    if not library.is_file() or ',' in str(library):
        raise ValueError('BLAS replacement must be a file with a linker-safe path')
    result = list(command)
    selected = []
    for index, argument in enumerate(result):
        name = Path(argument).name.lower()
        is_library = not argument.startswith('-') and name.startswith('libopenblas') and (
            name.endswith('.a') or '.so' in name or name.endswith('.dylib'))
        if argument == '-lopenblas' or is_library:
            selected.append(index)
    if not selected:
        raise RuntimeError('No explicit OpenBLAS link argument found; refusing an unreviewed fallback')
    for index in selected:
        result[index] = str(library)
    result.append('-Wl,-rpath,' + str(library.parent))
    return result


def main(argv=None):
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('mode', nargs='?')
    parser.add_argument('algebra_object', nargs='?', type=Path)
    parser.add_argument('--blas-library', type=Path)
    args = parser.parse_args(argv)
    build, output = args.build.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    flags = {}
    for line in (build / 'cpp/CMakeFiles/pea.dir/flags.make').read_text().splitlines():
        if ' = ' in line:
            key, value = line.split(' = ', 1)
            flags[key] = shlex.split(value)
    command = shlex.split((build / 'cpp/CMakeFiles/pea.dir/link.txt').read_text())
    if args.blas_library is not None:
        command = compatible_blas_link(command, args.blas_library)
    subprocess.run([
        command[0], *flags['CXX_DEFINES'], *flags['CXX_INCLUDES'], *flags['CXX_FLAGS'],
        *profile_definitions(args.mode), '-O2', '-fopenmp', '-std=c++20', '-c',
        str(root / 'src/tests/unit/r51_float_parallel_benchmark.cpp'),
        '-o', str(output / 'driver.o')], check=True)
    subprocess.run(['objcopy', '--redefine-sym', 'main=pea_original_main',
        str(build / 'cpp/CMakeFiles/pea.dir/pea/main.cpp.o'),
        str(output / 'original_main.o')], check=True)
    command[command.index('CMakeFiles/pea.dir/pea/main.cpp.o')] = str(output / 'original_main.o')
    if args.algebra_object is not None:
        # A real frozen production object, not the new core with an option off.
        command[command.index('CMakeFiles/pea.dir/common/algebra.cpp.o')] = str(args.algebra_object.resolve())
    command.insert(1, str(output / 'driver.o'))
    command[command.index('-o') + 1] = str(output / 'float_parallel_benchmark')
    subprocess.run(command, cwd=build / 'cpp', check=True)
    print(output / 'float_parallel_benchmark')


if __name__ == '__main__':
    main()
