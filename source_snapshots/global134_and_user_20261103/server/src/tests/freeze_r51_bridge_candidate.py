"""Freeze source/binary identity after the final successful build (not a test gate)."""
from pathlib import Path
import hashlib
import json
import shutil
import sys
import time

repo=Path(__file__).resolve().parents[2]
destination=Path(sys.argv[1]).resolve()
build=Path(sys.argv[2]).resolve()
assert destination.parent==Path('/home/rx/GINAN/bridge_continuity_20261020')
destination.mkdir(exist_ok=False)
binary=repo/'bin/pea'
assert binary.is_file()
shutil.copy2(binary,destination/'pea')
assert hashlib.sha256((destination/'pea').read_bytes()).digest()==hashlib.sha256(binary.read_bytes()).digest()
driver=repo/'bin/zhang_r51_bridge_event_replay_tests'
assert driver.is_file()
shutil.copy2(driver,destination/driver.name)
checkpoint_driver=destination.parent/'structural_checkpoint_verifier/structural_checkpoint_roundtrip'
assert checkpoint_driver.is_file()
shutil.copy2(checkpoint_driver,destination/checkpoint_driver.name)
files=[]
for path in sorted((repo/'src').rglob('*')):
    if path.is_file() and (path.suffix in ('.cpp','.hpp','.h','.py','.cmake') or path.name=='CMakeLists.txt'):
        relative=path.relative_to(repo)
        target=destination/'source'/relative
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(path,target)
        copied_hash=hashlib.sha256(target.read_bytes()).hexdigest()
        assert copied_hash==hashlib.sha256(path.read_bytes()).hexdigest(),'Source changed during freeze: '+str(relative)
        files.append(dict(path=str(relative),sha256=copied_hash,bytes=target.stat().st_size))
tests=destination/'tests'
tests.mkdir()
for source,name in ((build/'bridge_ctest.xml','bridge_ctest.xml'),
    (destination.parent/'cycle_slip_unchanged_regression.log','cycle_slip_unchanged_regression.log'),
    (destination.parent/'build_all_tests_guarded.log','build_all_tests_guarded.log'),
    (destination.parent/'ctest.log','ctest.log')):
    shutil.copy2(source,tests/name)
import subprocess
(tests/'configured_tests.json').write_bytes(subprocess.check_output(['ctest','--test-dir',str(build),'--show-only=json-v1']))
manifest=dict(created_unix=time.time(),repo=str(repo),scope='source_and_binary_freeze_not_scientific_acceptance',
    binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
    exact_replay_driver_sha256=hashlib.sha256((destination/driver.name).read_bytes()).hexdigest(),source_files=files,
    checkpoint_replay_driver_sha256=hashlib.sha256((destination/checkpoint_driver.name).read_bytes()).hexdigest(),
    test_evidence_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in tests.iterdir()},
    ratio=2,joint_nis=False,epoch_only=True,ar_started=False,
    required_before_formal_ar=['compiled_final_sources','all_required_tests_passed',
        'controlled_0206_event_verified','controlled_0237_event_verified',
        'real_exact_old_new_replays','checkpoint_restore_identity_verified'])
(destination/'manifest.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps(dict(binary_sha256=manifest['binary_sha256'],source_files=len(files),destination=str(destination))),flush=True)
