"""Link only a test main against already-finished PEA objects; never run PEA.

The build owner must finish the production build first. Output is fresh and
preserved on failure. No producer source/object is compiled or modified here.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('build', type=Path)
parser.add_argument('--output', required=True, type=Path)
parser.add_argument('--blas-library', type=Path)
parser.add_argument('--syntax-only', action='store_true')
parser.add_argument('--production-source-root', type=Path,
                    help='Explicit immutable source tree that owns the linked production objects; '
                         'CMake home must match this tree exactly. The test still comes from this runner tree.')
parser.add_argument('--legacy-objects-owner-regression-only', action='store_true',
                    help='RED-test only: omit the new dispatch API unavailable in old objects. '
                         'Such a run is never an acceptance-gate PASS.')
parser.add_argument('--legacy-r51-refresh-instrumented-red-only', action='store_true',
                    help='RED only: copy the old ppp_ambres object, expose three private symbols, '
                         'and replay its exact frozen old refresh loop. Original objects stay untouched.')
parser.add_argument('--legacy-gec-size-red-only', action='store_true',
                    help='RED only: use the current pure size-check header with unchanged V6 objects; '
                         'the old native importer must reject valid JSON over 256 MiB.')
args = parser.parse_args()
repo = Path(__file__).resolve().parents[2]
production_repo = (args.production_source_root or repo).resolve(strict=True)
build = args.build.resolve(strict=True)
output = args.output.resolve()
output.mkdir(parents=True, exist_ok=False)
cache = (build / 'CMakeCache.txt').read_text()
if ('CMAKE_HOME_DIRECTORY:INTERNAL=' + str(production_repo / 'src')) not in cache.splitlines():
    raise RuntimeError('Production objects belong to another source tree')
settings = {}
for line in (build / 'cpp/CMakeFiles/pea.dir/flags.make').read_text().splitlines():
    key, separator, value = line.partition(' = ')
    if separator:
        settings[key] = shlex.split(value)
link = shlex.split((build / 'cpp/CMakeFiles/pea.dir/link.txt').read_text())
if args.blas_library:
    library = args.blas_library.resolve(strict=True)
    selected = [i for i, item in enumerate(link) if 'openblas' in Path(item).name.lower()]
    if not selected:
        raise RuntimeError('No explicit OpenBLAS link; no fallback is permitted')
    for index in selected:
        link[index] = str(library)
    link.append('-Wl,-rpath,' + str(library.parent))
driver = repo / 'src/tests/unit/gec_ar_runtime_checkpoint.cpp'
source_paths = {'test/' + str(driver.relative_to(repo)): driver,
                'test/' + str(Path(__file__).resolve().relative_to(repo)): Path(__file__).resolve()}
limit_header = 'src/cpp/common/zhangGecArCheckpointLimits.hpp'
source_paths['test/support/' + limit_header] = repo / limit_header
if (production_repo / limit_header).is_file():
    source_paths['production/' + limit_header] = production_repo / limit_header
elif not args.legacy_gec_size_red_only:
    raise RuntimeError('Production shared GEC capacity header missing outside explicit old-size RED test')
for relative in [
    'src/cpp/pea/ppp_ambres_checkpoint.inc',
    'src/cpp/pea/ppp_ambres.cpp',
    'src/cpp/pea/zhangPppAr.hpp',
    'src/cpp/pea/zhangPppAr.cpp',
    'src/cpp/pea/outputs.cpp',
    'src/cpp/common/zhangCurrentEpochEvidence.hpp',
    'src/cpp/common/zhangMultiGnssInteger.hpp',
    'src/cpp/common/zhangMultiGnssFloatCheckpoint.hpp',
    'src/cpp/common/zhangIntegerDecisionProof.hpp',
]:
    source_paths['production/' + relative] = production_repo / relative

def source_hashes():
    return {label: hashlib.sha256(path.read_bytes()).hexdigest()
            for label, path in source_paths.items()}

source_before = source_hashes()
prefix = [link[0], *settings['CXX_DEFINES'], *settings['CXX_INCLUDES'],
          *settings['CXX_FLAGS'], '-O0', '-std=c++20', '-I' + str(output)]
if args.legacy_objects_owner_regression_only:
    prefix.append('-DZHANG_TEST_LEGACY_OBJECTS_NO_DISPATCH=1')
# Schema declarations only, for independent field contamination of the legacy
# wire. No validator, exporter, importer, or state-transition code is copied.
source = (production_repo / 'src/cpp/pea/ppp_ambres_checkpoint.inc').read_text().split('/** Multi-GNSS FLOAT checkpoints')[0]
def declaration(name, text=None):
    text = source if text is None else text
    match = re.search(r'^struct ' + re.escape(name) + r'\s*\{', text, re.MULTILINE)
    if not match:
        raise RuntimeError('Missing exact frozen struct declaration: ' + name)
    depth = 0
    for token in re.finditer(r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"|\{|\}', text[match.end()-1:]):
        if token.group() == '{': depth += 1
        elif token.group() == '}':
            depth -= 1
            if depth == 0:
                end = match.end()-1+token.end()
                if text[end] != ';': raise RuntimeError('Malformed schema declaration')
                return text[match.start():end+1]
    raise RuntimeError('Unterminated schema declaration')
names = re.findall(r'^struct (\w+)\s*\{', source, re.MULTILINE)
mirror = '#pragma once\nnamespace gec_test_ambres_wire {\n' + '\n\n'.join(declaration(name) for name in names) + '\n}\n'
(output / 'gec_private_ambres_wire_schema.hpp').write_text(mirror)
if args.syntax_only:
    raise SystemExit(subprocess.call([*prefix, '-fsyntax-only', str(driver)]))
objects = [build / 'cpp' / path for path in link if path.endswith('.o')]
def stamps():
    return {str(path.relative_to(build)): [path.stat().st_size, path.stat().st_mtime_ns]
            for path in objects}
before = stamps()
driver_object = output / 'gec_ar_runtime_checkpoint.o'
renamed_main = output / 'original_pea_main.o'
executable = output / 'gec_ar_runtime_checkpoint'
instrumentation = None
with (output / 'build.log').open('w') as log:
    subprocess.run([*prefix, '-c', str(driver), '-o', str(driver_object)],
                   stdout=log, stderr=subprocess.STDOUT, check=True)
    subprocess.run(['objcopy', '--redefine-sym', 'main=pea_original_main',
                    str(build / 'cpp/CMakeFiles/pea.dir/pea/main.cpp.o'), str(renamed_main)],
                   stdout=log, stderr=subprocess.STDOUT, check=True)
    link[link.index('CMakeFiles/pea.dir/pea/main.cpp.o')] = str(renamed_main)
    if args.legacy_r51_refresh_instrumented_red_only:
        # The original V5 loop is private and has no callable refresh API.
        # Expose symbols in a COPY, not the immutable producer object. The
        # bridge replays the exact extracted producer loop and calls the old
        # native store/exporter. No validator or integer solver is copied.
        cpp_path = production_repo / 'src/cpp/pea/ppp_ambres.cpp'
        cpp_source = cpp_path.read_text()
        loop_pattern = (r'if\(zhangR51Enabled\(\)\)for\(const auto& \[system,options\]:'
            r'acsConfig\.zhangFullRank\.sysOpts\)\s*'
            r'r51PhysicalArchives\[\{zhangAmbresRuntimeId\(kfState\),system\}\]='
            r'r51PhysicalStore\(kfState,system\);')
        matches = list(re.finditer(loop_pattern, cpp_source))
        if len(matches) != 1:
            raise RuntimeError('Expected exactly one frozen legacy terminal refresh loop')
        old_loop = matches[0].group()
        layout = declaration('ZhangR51PhysicalStore', cpp_source)
        preamble = '\n'.join(re.findall(r'^#include[^\n]*', cpp_source[:cpp_source.index('struct ZhangR51PhysicalStore')], re.MULTILINE))
        # cpp_source may contain includes only at file scope before this type;
        # checkpoint .inc inclusion occurs later and must never be duplicated.
        if 'ppp_ambres_checkpoint.inc' in preamble:
            raise RuntimeError('Unexpected checkpoint implementation in bridge preamble')
        symbols = ['_ZL16r51PhysicalStoreRK7KFState5E_Sys', '_ZL19r51PhysicalArchives',
            '_ZN12_GLOBAL__N_120zhangAmbresRuntimeIdERK7KFState']
        original_ambres = build / 'cpp/CMakeFiles/pea.dir/pea/ppp_ambres.cpp.o'
        copied_ambres = output / 'instrumented_copy_ppp_ambres.o'
        original_sha = hashlib.sha256(original_ambres.read_bytes()).hexdigest()
        subprocess.run(['objcopy', *['--globalize-symbol=' + symbol for symbol in symbols],
            str(original_ambres), str(copied_ambres)], stdout=log, stderr=subprocess.STDOUT, check=True)
        bridge = (preamble + '\nusing namespace std;\n' + layout + '\n' +
            'extern ZhangR51PhysicalStore r51PhysicalStore(const KFState&,E_Sys) asm("' + symbols[0] + '");\n' +
            'extern map<pair<string,E_Sys>,ZhangR51PhysicalStore> r51PhysicalArchives asm("' + symbols[1] + '");\n' +
            'extern string zhangAmbresRuntimeId(const KFState&) asm("' + symbols[2] + '");\n' +
            'void refreshZhangR51CurrentPhysicalArchives(const KFState& kfState) {\n' + old_loop + '\n}\n')
        bridge_path = output / 'legacy_r51_refresh_bridge.cpp'
        bridge_path.write_text(bridge)
        bridge_object = output / 'legacy_r51_refresh_bridge.o'
        subprocess.run([*prefix, '-c', str(bridge_path), '-o', str(bridge_object)],
            stdout=log, stderr=subprocess.STDOUT, check=True)
        link[link.index('CMakeFiles/pea.dir/pea/ppp_ambres.cpp.o')] = str(copied_ambres)
        link.insert(1, str(bridge_object))
        instrumentation = {'instrumented_RED_not_acceptance': True,
            'original_object_path': str(original_ambres), 'original_object_sha256': original_sha,
            'copied_object_path': str(copied_ambres),
            'copied_object_sha256': hashlib.sha256(copied_ambres.read_bytes()).hexdigest(),
            'globalized_symbols': symbols, 'frozen_source_path': str(cpp_path),
            'frozen_source_sha256': hashlib.sha256(cpp_path.read_bytes()).hexdigest(),
            'exact_old_loop': old_loop, 'exact_old_loop_sha256': hashlib.sha256(old_loop.encode()).hexdigest(),
            'layout_sha256': hashlib.sha256(layout.encode()).hexdigest(),
            'bridge_path': str(bridge_path), 'bridge_sha256': hashlib.sha256(bridge_path.read_bytes()).hexdigest(),
            'bridge_object_sha256': hashlib.sha256(bridge_object.read_bytes()).hexdigest()}
    link.insert(1, str(driver_object))
    link[link.index('-o') + 1] = str(executable)
    subprocess.run(link, cwd=build / 'cpp', stdout=log, stderr=subprocess.STDOUT, check=True)
if before != stamps() or source_before != source_hashes():
    raise RuntimeError('Production objects or contract sources changed during link; discard this test result')
env = {key: value for key, value in os.environ.items() if not key.startswith('ZHANG_')}
env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
with (output / 'runtime.log').open('w') as log:
    result = subprocess.run([str(executable)], cwd=output, env=env,
                            stdout=log, stderr=subprocess.STDOUT, timeout=120)
marker = 'GEC_AR_ARCHIVE_RESULT '
records = [json.loads(line.partition(marker)[2]) for line in (output / 'runtime.log').read_text().splitlines()
           if line.startswith(marker)]
report = {'native_exit_code': result.returncode,
          'production_source_root': str(production_repo),
          'test_source_root': str(repo),
          'production_cmake_home_exact': True,
          'legacy_objects_owner_regression_only': args.legacy_objects_owner_regression_only,
          'legacy_r51_refresh_instrumented_red_only': args.legacy_r51_refresh_instrumented_red_only,
          'legacy_gec_size_red_only': args.legacy_gec_size_red_only,
          'RED_instrumentation': instrumentation,
          'driver_sha256': hashlib.sha256(driver.read_bytes()).hexdigest(),
          'runner_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'executable_sha256': hashlib.sha256(executable.read_bytes()).hexdigest(),
          'production_object_stamps_unchanged': before == stamps(),
          'production_object_stamps': before,
          'contract_source_sha256': source_before,
          'contract_sources_unchanged': source_before == source_hashes(),
          'legacy_wire_schema_sha256': hashlib.sha256(
              (output / 'gec_private_ambres_wire_schema.hpp').read_bytes()).hexdigest(),
          'result': records[-1] if records else None,
          'passed': result.returncode == 0 and len(records) == 1 and
              not args.legacy_objects_owner_regression_only and not args.legacy_r51_refresh_instrumented_red_only and
              not args.legacy_gec_size_red_only and
              records[0]['status'] == 'PASS' and before == stamps() and
              source_before == source_hashes()}
(output / 'result.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
raise SystemExit(0 if report['passed'] else 1)
