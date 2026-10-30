#!/usr/bin/env python3
"""Freeze the complete source and a private producer identity; not a test pass."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[2]
    binary = args.binary.resolve(strict=True)
    destination = args.output.resolve()
    destination.mkdir(parents=True, exist_ok=False)
    files = sorted(path for base in ('src', 'experiments', 'cmake')
                   for path in (repo / base).rglob('*') if path.is_file()
                   and '__pycache__' not in path.parts and path.suffix not in ('.pyc', '.o', '.a', '.so'))
    for name in ('CMakeLists.txt', 'CMakePresets.json'):
        if (repo / name).is_file(): files.append(repo / name)
    identities = []
    archive = destination / 'source.zip'
    with zipfile.ZipFile(archive, 'x', compression=zipfile.ZIP_DEFLATED) as output:
        for path in files:
            identity = {'path': path.relative_to(repo).as_posix(),
                        'bytes': path.stat().st_size, 'sha256': sha(path)}
            output.write(path, identity['path'])
            if identity['sha256'] != sha(path):
                raise RuntimeError('Source changed during freeze: ' + identity['path'])
            identities.append(identity)
    with zipfile.ZipFile(archive) as saved:
        for row in identities:
            if hashlib.sha256(saved.read(row['path'])).hexdigest() != row['sha256']:
                raise RuntimeError('Archive content differs: ' + row['path'])
    manifest = {'schema': 1, 'binary': str(binary), 'binary_sha256': sha(binary),
                'archive_sha256': sha(archive), 'source_files': identities,
                'scope': 'Immutable source/producer binding only; not test or scientific acceptance'}
    (destination / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps({key: value for key, value in manifest.items() if key != 'source_files'}))
    print('source_files=' + str(len(identities)))


if __name__ == '__main__':
    main()
