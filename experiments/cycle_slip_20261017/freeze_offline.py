"""Freeze a qualified offline experiment. Does not launch PEA or modify inputs."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile


def sha(path):
    value=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):value.update(block)
    return value.hexdigest()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--validation',type=Path,required=True)
    parser.add_argument('--frozen',type=Path,required=True)
    parser.add_argument('--run-root',type=Path,required=True)
    parser.add_argument('--restart-after-passed-tests',action='store_true')
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[2]
    report=json.loads(args.validation.read_text())
    if not args.restart_after_passed_tests:raise SystemExit('The user must have requested restart after passed tests')
    if not report['all_execution_checks_passed']:raise SystemExit('Validation failed; no freeze/launch')
    if sha(root/'bin/pea')!=report['binary_sha256']:raise SystemExit('PEA changed after validation')
    for name,digest in report['source_hashes'].items():
        if sha(root/name)!=digest:raise SystemExit('Source changed after validation: '+name)
    for path,digest in report['input_hashes'].items():
        if sha(Path(path))!=digest:raise SystemExit('Original input changed after validation: '+path)
    if args.run_root.exists():raise SystemExit('Run output already exists')
    if args.run_root.parent.resolve()!=Path('/home/rx/GINAN/local_timing') or not args.run_root.name.startswith('r51_20261001_osb_branch_repair_ar1h_'):
        raise SystemExit('Unexpected run output scope')
    evidence=args.validation.parent
    injection=json.loads((evidence/'injected_independent_holdout.log').read_text().splitlines()[-1])
    null=json.loads((evidence/'null_independent_holdout.log').read_text().splitlines()[-1])
    if injection['detected_cases']!=injection['injected_cases'] or injection['wrong_integer_certificates'] or null['retirements']:
        raise SystemExit('Scientific readiness counters failed')
    args.frozen.mkdir(parents=True,exist_ok=False)
    shutil.copy2(root/'bin/pea',args.frozen/'pea')
    shutil.copy2(Path(__file__).with_name('launch_offline.sh'),args.frozen/'launch_offline.sh')
    shutil.copy2(args.validation,args.frozen/'validation_result.json')
    original=Path('/home/rx/GINAN/frozen-r51-satellite-fix-20260918/config')
    shutil.copytree(original,args.frozen/'config')
    case=args.frozen/'config/case.yaml'
    text=case.read_text()
    if text.count('\noutputs:\n')!=1:raise SystemExit('Unexpected case.yaml layout')
    case.write_text(text.replace('\noutputs:\n',
        '\noutputs:\n  trace:\n    level: 2\n    output_receivers: true\n',1))
    with tarfile.open(args.frozen/'source.tar.gz','w:gz') as archive:
        archive.add(root/'src',arcname='src')
        archive.add(root/'experiments/cycle_slip_20261017',arcname='experiments/cycle_slip_20261017')
    checkpoint=Path('/mnt/d/GINAN_R20/zhang_r51_ratio3_float2h_ar1h_2024199_180_20260917_checkpoints_source_before_ratioonly/E29-20240717015930-e240-rb06f1454a2e1/checkpoint.bundle')
    if sha(checkpoint)!=report['checkpoint_sha256']:raise SystemExit('Original checkpoint changed')
    launcher=Path('/home/rx/GINAN/frozen-osb-ratio2-domain-budget-fix-20261001-v5/launcher.sh')
    manifest={
        'status':'FROZEN_NOT_STARTED','binary_sha256':sha(args.frozen/'pea'),
        'source_archive_sha256':sha(args.frozen/'source.tar.gz'),
        'validation_report':str(args.validation),'run_root':str(args.run_root),
        'base_launcher':str(launcher),'base_launcher_sha256':sha(launcher),
        'checkpoint':str(checkpoint),'checkpoint_sha256':report['checkpoint_sha256'],
        'data_span_gpst':['2024-07-17 02:00:00','2024-07-17 02:59:30'],
        'ratio':2,'epoch_only_ar':True,'network_joint_nis':'SKIPPED_RATIO_ONLY',
        'slip_qc':'MATCHED_PEER_WINDOW','offline_raw_qc_warmup_seconds':1920,
        'qualified_gf_assistance':'PAST_ONLY_SINGLE_AND_FROZEN_EVENT_WINDOW',
        'early_auxiliary_policy':'QUARANTINE_ONLY_NOT_INTEGER_CERTIFICATE',
        'early_stage_trace':{'0':'NONE','1':'CURRENT_QUALIFIED_GF','2':'FROZEN_EVENT_GF','3':'FIXED_LAG_DOPPLER'},
        'instantaneous_test_family':{'without_doppler':2,'with_doppler':6},
        'detector_history_bounds':{'accepted':32,'screening':96},
        'long_tail_regression_gates':report['long_tail_regression_gates'],
        'input_hashes':report['input_hashes'],
        'restart_authorization':'USER_REQUEST_AFTER_TESTS_PASSED',
        'performance_gates':report['performance_gates'],
        'measured_injected_delay_max_seconds':injection['max_detection_delay_s'],
        'measured_injected_delay_p95_seconds':injection['p95_detection_delay_s'],
        'injected_test_summary':injection,'null_test_summary':null,
        'config_hashes':{p.name:sha(p) for p in sorted((args.frozen/'config').glob('*.yaml'))},
        'trace_change':'Receiver TRACE enabled, level 2, for causal audit; timing not strictly comparable to V5 level 1',
        'limitations':report['limitations']+[
            'An injected slip may contaminate the float model before delayed quarantine',
            'Detection/quarantine is not successful integer identification or repair',
            'Local 2-4 dimensional covariance-weighted slip QC remains enabled; network joint NIS remains disabled',
            'Real-time safety and product continuity are not established by these QC tests'],
        'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
        'git_status':subprocess.check_output(['git','status','--porcelain'],cwd=root,text=True)}
    (args.frozen/'manifest.json').write_text(json.dumps(manifest,indent=2))
    print(json.dumps({'frozen':str(args.frozen),'run_root':str(args.run_root),
        'status':manifest['status'],'max_detection_delay_s':manifest['measured_injected_delay_max_seconds']}))


if __name__=='__main__':main()
