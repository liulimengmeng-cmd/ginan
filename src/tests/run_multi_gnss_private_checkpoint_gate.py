"""Compile only a test driver, then link stable production PEA objects.

Never recompiles production source or starts a GNSS solver. The parent must
finish its build before invocation. --syntax-only is independent of objects.
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
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--syntax-only', action='store_true')
parser.add_argument('--producer-syntax-only', nargs='?', const='src/cpp/pea/zhangPppAr.cpp')
parser.add_argument('--blas-library', type=Path)
args = parser.parse_args()
repo = Path(__file__).resolve().parents[2]
build = args.build.resolve()
output = args.output.resolve()
output.mkdir(parents=True, exist_ok=True)
cache = (build / 'CMakeCache.txt').read_text()
assert 'CMAKE_HOME_DIRECTORY:INTERNAL=' + str(repo / 'src') in cache
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
driver = repo / 'src/tests/unit/multi_gnss_private_checkpoint_gate.cpp'
prefix = [link[0], *settings['CXX_DEFINES'], *settings['CXX_INCLUDES'],
          *settings['CXX_FLAGS'], '-O0', '-std=c++20', '-I' + str(output)]

# A test-only schema mirror permits field-level wire contamination without
# exposing a producer API or compiling any producer implementation in the
# driver. This mechanically copies DTO declarations and serialization only,
# never the gate, export/import, state mutation, or validation algorithms.
# The assertions below are independent expected-policy checks, not a second
# invocation of a copied reference implementation.
producer_source = repo / 'src/cpp/pea/zhangPppAr.cpp'
source = producer_source.read_text()
names = ['ProductKey', 'E27RawNoiseRow', 'PromotionEvidenceKey', 'PromotionEvidence',
         'ProductHistoryKey', 'UserReferenceKey', 'UserReferenceState', 'UserDualReferenceKey',
         'ZhangCapturedFactorEventCheckpoint', 'ZhangCapturedPhysicalTargetCheckpoint',
         'ZhangCapturedUnresolvedIntegerDatumCheckpoint', 'ZhangCapturedRetainedTargetBlockCheckpoint',
         'ZhangFactorCaptureCheckpoint', 'ZhangPhaseContinuityCheckpoint', 'GlobalContinuityCheckpoint',
         'ZhangPendingSnapshotPinsCheckpoint', 'E27NoiseSensitivityCheckpoint',
         'E27JointNoiseRuntimeCheckpoint', 'ProductHistoryCheckpoint',
         'ProductPublicationIdentityCheckpoint', 'ZhangGraphEdgeCheckpoint',
         'ZhangProductIntegerFunctionalCheckpoint', 'ZhangProductIntegerTransitionCheckpoint',
         'ZhangFullJointStructuralChainCheckpoint', 'ZhangPendingProductTransitionCheckpoint',
         'E27RawNoiseRowCheckpoint', 'ZhangHybridRealGaugeCheckpointEntry', 'ZhangPppArCheckpointEnvelope']

def declaration(name, text=source):
    match = re.search(r'^struct ' + re.escape(name) + r'\s*\{', text, re.MULTILINE)
    if not match:
        raise RuntimeError('Missing controlled test DTO: ' + name)
    # Tokenization ignores comments and string contents while balancing C++
    # scopes. No arbitrary byte search/patching of the checkpoint is used.
    tokens = re.finditer(r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"|\{|\}',
                         text[match.end() - 1:])
    depth = 0
    for token in tokens:
        if token.group() == '{':
            depth += 1
        elif token.group() == '}':
            depth -= 1
            if depth == 0:
                end = match.end() - 1 + token.end()
                assert text[end] == ';'
                return text[match.start():end + 1]
    raise RuntimeError('Unterminated controlled DTO: ' + name)

serializers = source[source.index('namespace boost::serialization'):source.index('} // namespace boost::serialization')]
mirror = '''// Generated schema declarations only; not a reference gate.\n#pragma once
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/deque.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/utility.hpp>
#include "common/zhangFactorCapture.hpp"
#include "common/zhangPhaseContinuity.hpp"
#include "common/zhangHybridService.hpp"
using std::map; using std::set; using std::string; using std::vector;
'''
mirror += serializers + '} // namespace boost::serialization\nnamespace checkpoint_test_wire {\n'
mirror += '\n\n'.join(declaration(name) for name in names) + '\n}\n'
(output / 'private_checkpoint_wire_schema.hpp').write_text(mirror)
ambres_source = repo / 'src/cpp/pea/ppp_ambres_checkpoint.inc'
ambres_text = ambres_source.read_text().split('/** Multi-GNSS FLOAT checkpoints')[0]
ambres_names = re.findall(r'^struct (\w+)\s*\{', ambres_text, re.MULTILINE)
ambres_mirror = '#pragma once\nnamespace checkpoint_test_ambres_wire {\n'
ambres_mirror += '\n\n'.join(declaration(name, ambres_text) for name in ambres_names) + '\n}\n'
(output / 'private_ambres_wire_schema.hpp').write_text(ambres_mirror)
graph_source = repo / 'src/cpp/pea/zhangReference.cpp'
graph_text = graph_source.read_text()
graph_names = ['ZhangGraphCheckpointSatellite', 'ZhangGraphCheckpointEdge',
               'ZhangGraphCheckpointBasis', 'ZhangGraphCheckpointEdgeHistory',
               'ZhangGraphCheckpointArcObservationHistory', 'ZhangGraphCheckpointOutageState',
               'ZhangGraphCheckpointRuntimeState', 'ZhangGraphCheckpointRuntimeV4',
               'ZhangGraphCheckpointPayloadV4']
graph_magic = re.search(r'constexpr char ZHANG_GRAPH_CHECKPOINT_MAGIC\[\]\s*=\s*"[^"]+";', graph_text).group()
graph_mirror = '#pragma once\n#include "pea/zhangReference.hpp"\nnamespace checkpoint_test_graph_wire {\n' + graph_magic + '\n'
graph_mirror += '\n\n'.join(declaration(name, graph_text) for name in graph_names) + '\n}\n'
(output / 'private_graph_wire_schema.hpp').write_text(graph_mirror)
(output / 'wire_schema_provenance.json').write_text(json.dumps({
    'source_sha256': hashlib.sha256(producer_source.read_bytes()).hexdigest(),
    'declarations': names, 'copied_core_functions': False,
    'ambres_source_sha256': hashlib.sha256(ambres_source.read_bytes()).hexdigest(),
    'ambres_declarations': ambres_names,
    'graph_source_sha256': hashlib.sha256(graph_source.read_bytes()).hexdigest(),
    'graph_declarations': graph_names,
    'purpose': 'Wire field injection only; policy assertions use actual production adapters'
}, indent=2) + '\n')
if args.producer_syntax_only:
    assert args.producer_syntax_only in {'src/cpp/pea/zhangPppAr.cpp', 'src/cpp/pea/ppp_ambres.cpp'}
    raise SystemExit(subprocess.call([*prefix, '-fsyntax-only',
        str(repo / args.producer_syntax_only)]))
if args.syntax_only:
    raise SystemExit(subprocess.call([*prefix, '-fsyntax-only', str(driver)]))

objects = [build / 'cpp' / item for item in link if item.endswith('.o')]
def stamps():
    return {str(path.relative_to(build)): [path.stat().st_size, path.stat().st_mtime_ns]
            for path in objects}
before = stamps()
driver_object = output / 'private_checkpoint_gate_driver.o'
renamed_main = output / 'original_pea_main.o'
executable = output / 'private_checkpoint_gate_contract'
with (output / 'build.log').open('w') as log:
    subprocess.run([*prefix, '-c', str(driver), '-o', str(driver_object)],
                   stdout=log, stderr=subprocess.STDOUT, check=True)
    subprocess.run(['objcopy', '--redefine-sym', 'main=pea_original_main',
                    str(build / 'cpp/CMakeFiles/pea.dir/pea/main.cpp.o'), str(renamed_main)],
                   stdout=log, stderr=subprocess.STDOUT, check=True)
    link[link.index('CMakeFiles/pea.dir/pea/main.cpp.o')] = str(renamed_main)
    link.insert(1, str(driver_object))
    link[link.index('-o') + 1] = str(executable)
    subprocess.run(link, cwd=build / 'cpp', stdout=log, stderr=subprocess.STDOUT, check=True)
assert before == stamps(), 'Production objects changed during test link; discard this result'
env = {key: value for key, value in os.environ.items() if not key.startswith('ZHANG_')}
env.update(OMP_NUM_THREADS='8', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
with (output / 'runtime.log').open('w') as log:
    result = subprocess.run([str(executable)], cwd=output, env=env,
                            stdout=log, stderr=subprocess.STDOUT, timeout=90)
runtime_text = (output / 'runtime.log').read_text()
marker = 'PRIVATE_CHECKPOINT_GATE_RESULT '
records = [json.loads(line.partition(marker)[2]) for line in runtime_text.splitlines() if marker in line]
report = {'driver_sha256': hashlib.sha256(driver.read_bytes()).hexdigest(),
          'executable_sha256': hashlib.sha256(executable.read_bytes()).hexdigest(),
          'exit_code': result.returncode, 'production_object_stamps_unchanged': before == stamps(),
          'selected_blas_library': str(args.blas_library.resolve()) if args.blas_library else 'original-production-CMake-link',
          'production_object_stamps': before, 'result': records[-1] if records else None,
          'passed': result.returncode == 0 and bool(records) and records[-1]['status'] == 'PASS'}
(output / 'result.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({key: report[key] for key in ('passed', 'exit_code', 'result')}, indent=2))
raise SystemExit(0 if report['passed'] else 1)
