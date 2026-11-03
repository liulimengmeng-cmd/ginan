"""Link-command fixtures only; no C++ compiler/solver is launched."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import run_r51_float_parallel as builder


class GenericBlasBuildContracts(unittest.TestCase):
    def test_original_baseline_and_frozen_profile_definitions(self):
        self.assertEqual(builder.profile_definitions(None), ['-DFLOAT_THREAD_PROFILE=1'])
        self.assertEqual(builder.profile_definitions('baseline'), [])
        self.assertEqual(builder.profile_definitions('frozen'),
            ['-DFLOAT_THREAD_PROFILE=1', '-DFLOAT_FROZEN_PROFILE=1'])

    def test_all_static_and_dynamic_libraries_replaced_with_explicit_rpath(self):
        with tempfile.TemporaryDirectory() as temp:
            library = Path(temp) / 'libopenblasp-r0.3.20.so'
            library.write_bytes(b'fixture-no-native-code')
            command = ['/usr/bin/c++', 'one.o', '/old/libopenblas.a', '/old/libopenblas.so.0',
                       '-Wl,-rpath,/do-not-change/openblas', '-o', 'pea']
            original = list(command)
            actual = builder.compatible_blas_link(command, library)
            self.assertEqual(actual[2:4], [str(library.resolve())] * 2)
            self.assertEqual(actual[-1], '-Wl,-rpath,' + str(library.parent.resolve()))
            self.assertIn('-Wl,-rpath,/do-not-change/openblas', actual)
            self.assertEqual(command, original)

    def test_absent_explicit_blas_does_not_silently_fall_back(self):
        with tempfile.TemporaryDirectory() as temp:
            library = Path(temp) / 'libgeneric.so'
            library.write_bytes(b'fixture')
            with self.assertRaisesRegex(RuntimeError, 'No explicit OpenBLAS'):
                builder.compatible_blas_link(['c++', 'openblas_helper.o', '-DGINAN_USE_OPENBLAS=1'], library)
            with self.assertRaises(FileNotFoundError):
                builder.compatible_blas_link(['c++', '/old/libopenblas.a'], library.parent / 'missing.so')

    def fixture_build(self, root):
        directory = root / 'build/cpp/CMakeFiles/pea.dir'
        directory.mkdir(parents=True)
        (directory / 'flags.make').write_text('CXX_DEFINES = -DTEST=1\nCXX_INCLUDES = -I/example\nCXX_FLAGS = -O2\n')
        (directory / 'link.txt').write_text('c++ CMakeFiles/pea.dir/pea/main.cpp.o '
            'CMakeFiles/pea.dir/common/algebra.cpp.o /old/libopenblas.a -o pea\n')
        return root / 'build'

    def test_original_positional_frozen_algebra_contract_and_new_option_combine(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            build = self.fixture_build(root)
            library = root / 'libopenblas.so'
            library.write_bytes(b'fixture')
            algebra = root / 'frozen_algebra.o'
            algebra.write_bytes(b'fixture')
            output = root / 'output'
            with patch.object(builder.subprocess, 'run') as run:
                builder.main([str(build), str(output), 'frozen', str(algebra), '--blas-library', str(library)])
            self.assertEqual(run.call_count, 3)
            compile_command = run.call_args_list[0].args[0]
            self.assertIn('-DFLOAT_THREAD_PROFILE=1', compile_command)
            self.assertIn('-DFLOAT_FROZEN_PROFILE=1', compile_command)
            link_command = run.call_args_list[2].args[0]
            self.assertIn(str(algebra.resolve()), link_command)
            self.assertIn(str(library.resolve()), link_command)
            self.assertIn('-Wl,-rpath,' + str(library.parent.resolve()), link_command)
            self.assertNotIn('/old/libopenblas.a', link_command)

    def test_no_new_option_keeps_original_default_and_baseline_linking(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            build = self.fixture_build(root)
            for mode in (None, 'baseline'):
                with patch.object(builder.subprocess, 'run') as run:
                    arguments = [str(build), str(root / ('default' if mode is None else mode))]
                    builder.main(arguments + ([] if mode is None else [mode]))
                compile_command = run.call_args_list[0].args[0]
                self.assertEqual('-DFLOAT_THREAD_PROFILE=1' in compile_command, mode is None)
                self.assertIn('/old/libopenblas.a', run.call_args_list[2].args[0])


if __name__ == '__main__':
    unittest.main()
