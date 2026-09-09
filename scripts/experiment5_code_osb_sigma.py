#!/usr/bin/env python3
"""Controlled daily CODE phase OSB prior sensitivity; never alters original inputs.

Only GPS satellite L1C/L2W sigma bytes (92:103) in derived experimental BIA
fixtures change. Means and all remaining bytes must be identical. These are
assumed independent satellite/signal priors, persistent and shared across sites,
NOT CODE formal uncertainties. The same frozen baseline PEA processes all runs.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time
from datetime import datetime, timezone

C = 299792458.0
EXPECTED_PEA = '3ab177d7ee7c06936d02384368fe0f090ccee0a16d4e6b6049fadde8b486b23c'
EXPECTED_BIA = 'fb99686fff296688920df31e92720d3f66b45c9991aa2ff42f74e732b1542483'
REPO = Path('/mnt/d/tec/ginan-main-code-products')
DOC = REPO / 'Docs/stecCovarianceFeasibility'
DATA = Path('/home/rx/GINAN/inputData')
PEA = Path('/mnt/d/tec/ginan-main-full-rank-ar/bin/pea')

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()

def write_json(path, value):
    with Path(path).open('x', encoding='utf-8') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')

def perturb(data, sigma_m):
    assert sigma_m > 0
    output, changed = [], []
    for i, line in enumerate(data.splitlines(keepends=True)):
        target = (line[1:5] == b'OSB ' and line[11:12] == b'G'
                  and line[15:19] == b'    ' and line[25:28] in (b'L1C', b'L2W'))
        if target:
            assert line[65:68] == b'ns ' and float(line[92:103]) == 0
            value = f'{sigma_m / C * 1e9:11.7f}'.encode()
            assert len(value) == 11
            new = line[:92] + value + line[103:]
            assert new[:92] == line[:92] and new[103:] == line[103:]
            changed.append({'line': i+1, 'satellite': line[11:14].decode(),
                            'signal': line[25:28].decode(), 'sigma_ns': float(value),
                            'effective_sigma_m': float(value)*C/1e9})
            output.append(new)
        else:
            output.append(line)
    assert len(changed) == 64
    return b''.join(output), changed

def prepare(root):
    # This linked worktree was created by Windows Git; its .git path is Windows-native.
    gitdir = (REPO/'.git').read_text().strip().removeprefix('gitdir: ')
    if len(gitdir) > 2 and gitdir[1] == ':':
        gitdir = '/mnt/' + gitdir[0].lower() + gitdir[2:]
    source_commit = subprocess.check_output(['git', '--git-dir', gitdir, 'rev-parse', 'HEAD'], text=True).strip()
    assert source_commit == '91c1ceae9ecf1192419d1db83841d1751ef7cb37'
    root.mkdir(parents=True, exist_ok=False)
    assert sha(PEA) == EXPECTED_PEA
    shutil.copy2(PEA, root / 'pea')
    manifest = json.loads((DOC/'experiment_1_2024_input_manifest.json').read_text())
    files = [f for f in manifest['files'] if 'WUM' not in f['relative_path'] and '.OBX' not in f['relative_path'].upper()]
    code = json.loads((DOC/'experiment_4_code_rapid_run_receipt.json').read_text())
    for name, digest in code['code_products_sha256'].items():
        files.append({'role': 'CODE rapid', 'relative_path': 'products/code_2024199/'+name, 'sha256': digest})
    for f in files:
        p = DATA/f['relative_path']
        assert sha(p) == f['sha256'], p
        f['size'] = p.stat().st_size
        f['absolute_path'] = str(p)
    write_json(root/'all_inputs.json', files)
    bia = DATA/'products/code_2024199/COD0OPSRAP_20241990000_01D_01D_OSB.BIA'
    assert sha(bia) == EXPECTED_BIA
    plans = []
    for mm in (3, 10):
        fixture = root/f'assumed_sigma_{mm}mm.BIA'
        content, changes = perturb(bia.read_bytes(), mm/1000)
        fixture.write_bytes(content)
        write_json(root/f'fixture_{mm}mm_receipt.json', {
            'original_path': str(bia), 'original_sha256': sha(bia),
            'fixture_path': str(fixture), 'fixture_sha256': sha(fixture),
            'assumed_sigma_m': mm/1000, 'changed_fields': changes,
            'only_bytes_92_103_of_64_GPS_phase_OSB_lines_changed': True,
            'published_CODE_product': False})
        # Standalone copy prevents GINAN's append-only input lists loading two BIAs.
        base = (DOC/'experiment_4_code_rapid_base.yaml').read_text()
        old = 'code_2024199/COD0OPSRAP_20241990000_01D_01D_OSB.BIA'
        assert base.count(old) == 1
        base_path = root/f'base_{mm}mm.yaml'
        base_path.write_text(base.replace(old, str(fixture)))
        for mode, original in [('ar', '4c_2880e_code_rapid_persistent_phase_osb'),
                               ('float', '4d_2880e_code_rapid_persistent_phase_osb_float')]:
            label = f'sigma_{mm}mm_{mode}'
            cfg = (DOC/f'experiment_{original}.yaml').read_text()
            cfg = cfg.replace('Docs/stecCovarianceFeasibility/experiment_4_code_rapid_base.yaml', str(base_path))
            cfg_path = root/f'{label}.yaml'
            cfg_path.write_text(cfg)
            out = root/label
            plans.append({'label': label, 'assumed_sigma_m': mm/1000, 'mode': mode,
                          'config': str(cfg_path), 'config_sha256': sha(cfg_path),
                          'base_sha256': sha(base_path), 'output': str(out),
                          'argv': [str(root/'pea'), '-y', str(cfg_path), '-n', '2880',
                                   '-a', f'EXP1_DATA_ROOT:{DATA}', '-a', f'EXP1_OUTPUT_ROOT:{out}']})
    write_json(root/'plan.json', {'created_utc': datetime.now(timezone.utc).isoformat(),
        'source_commit': source_commit,
        'binary_sha256': sha(root/'pea'), 'runs': plans,
        'sigma_grid_predeclared_m': [0.003, 0.010], 'zero_baseline': 'experiment 4c/4d',
        'covariance': 'Independent daily satellite/signal priors; shared across all receivers; Q=0.',
        'analysis_window': '2024-07-17 14:00:00 through 23:59:30',
        'claim': 'Sensitivity experiment, not calibrated uncertainty or certified fixing.'})
    print(root, flush=True)

def run(root, label, epochs):
    plan = json.loads((root/'plan.json').read_text())
    r = next(r for r in plan['runs'] if r['label'] == label)
    assert sha(root/'pea') == plan['binary_sha256']
    assert sha(r['config']) == r['config_sha256']
    suffix = '' if epochs == 2880 else f'_smoke_{epochs}'
    out = Path(r['output']+suffix)
    out.mkdir(exist_ok=False)
    argv = r['argv'].copy()
    argv[argv.index('-n')+1] = str(epochs)
    argv[-1] = f'EXP1_OUTPUT_ROOT:{out}'
    env = os.environ.copy()
    env.update(LD_LIBRARY_PATH='/home/rx/.local/boost-1.82/lib', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    started = datetime.now(timezone.utc).isoformat()
    t = time.monotonic()
    with (root/f'{label}{suffix}.stdout.log').open('x') as stdout, (root/f'{label}{suffix}.stderr.log').open('x') as stderr:
        p = subprocess.Popen(argv, cwd=REPO, env=env, stdout=stdout, stderr=stderr)
        write_json(root/f'{label}{suffix}.started.json', {'pid': p.pid, 'argv': argv, 'started_utc': started,
                   'environment_overrides': {k: env[k] for k in ['LD_LIBRARY_PATH','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']}})
        rc = p.wait()
    write_json(root/f'{label}{suffix}.finished.json', {'returncode': rc, 'elapsed_seconds': time.monotonic()-t,
                'finished_utc': datetime.now(timezone.utc).isoformat(), 'binary_sha256_after': sha(root/'pea')})
    print(label, epochs, rc, flush=True)
    if rc:
        raise SystemExit(rc)

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('action', choices=['prepare','run'])
    p.add_argument('root', type=Path)
    p.add_argument('--label')
    p.add_argument('--epochs', type=int, default=2880)
    a = p.parse_args()
    prepare(a.root) if a.action == 'prepare' else run(a.root, a.label, a.epochs)
