"""Prepare the exact original-e240 local datum-repair AR hour; never launch."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path

ART=Path('/mnt/c/Users/rx/Documents/GINAN/artifacts/checkpoint_authority_fix_20261104_01')
ROOT='/home/rx/GINAN/checkpoint_authority_fix_20261104_01'
BINARY=ROOT+'/candidate/pea'
AUDITOR=ROOT+'/candidate/us20_e240_nine_section_audit'
BASE=ART.parent/'ar_p0_p4_5h_20261102_01/harness/prepare_e240_5h_plan.py'
BASE_SHA='987b95b8b668c2eab4b38a0967682f60c1d10293f458c6ddd3f6908623073360'
POLICY='US20_GLOBAL134_E240_CODEC_FIX_V2'
SOURCE='/home/rx/GINAN/us20_global134_osb5h_20261031_v1/run/checkpoints/E29-20240717015930-e240-r32f626de9463'
SOURCE_ID='E29-20240717015930-e240-r32f626de9463'
SOURCE_BINARY='/home/rx/GINAN/us20_global134_osb5h_20261031_v1/candidate/pea'
SOURCE_BINARY_SHA='263db653de36e433d2458de7a34064ee3090fb6aa834324b5fb5e5e78ba5da16'
SOURCE_BUNDLE_SHA='abe8b436a814cd2f9dab87b979675c1328ba85aa5a1cc797cbdcc91a28189db5'
SOURCE_MANIFEST_SHA='690ca1d4559dfb312d3b30138081d8830fb9ed1c4506489d3cd540f86d34f517'
SOURCE_CONFIG_SHA='752fa1a4ef93d4f6529767d111b80bcfbe826207d698b063322f124feb23946c'
SOURCE_INPUT_SHA='e325b47f85b6061b0859bf9af7e9256c7f6bd88a1cf2830099cd4faca06f5da4'
FLAG='ZHANG_US20_E240_TRANSPORT_CONTINUATION'


def sha(raw):return hashlib.sha256(raw).hexdigest()


def host_path(path):
    p=str(path).replace('\\','/')
    return Path(p[5]+':/'+p[7:]) if os.name=='nt' and p.startswith('/mnt/') else Path(p)


def base_module():
    p=host_path(BASE)
    if sha(p.read_bytes())!=BASE_SHA:raise RuntimeError('Frozen canonical codec helper changed')
    s=importlib.util.spec_from_file_location('_codec_replay_canonical_base',p)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


P=base_module()
PROFILES={
    'formal':dict(run=ROOT+'/run',end='2024-07-17 02:59:30',max_epoch=360,
        captures=[(300,'2024-07-17 02:29:30'),(360,'2024-07-17 02:59:30')],
        overlay_name='100_checkpoint_authority_ar1h_overlay.yaml',
        description='Checkpoint authority repair direct AR1h',ar_start=1405216800,ar_first_epoch=241,actual_ar_epochs=120),
    'replay_e300':dict(run=ROOT+'/validation_e300/run',end='2024-07-17 02:29:30',max_epoch=300,
        captures=[(300,'2024-07-17 02:29:30')],
        overlay_name='100_checkpoint_authority_e300_replay_overlay.yaml',
        description='Checkpoint authority repair e300 replay',ar_start=1405216800,ar_first_epoch=241,actual_ar_epochs=60)}


def overlay(profile):
    p=PROFILES[profile];run=p['run']
    captures=''.join("        - '"+epoch+"'\n" for _,epoch in p['captures'])
    return ("processing_options:\n  epoch_control:\n"
        f"    end_epoch: '{p['end']}'\n    max_epochs: {p['max_epoch']}\n"
        "  gnss_general:\n    zhang_pppar:\n"
        f"      checkpoint_output_directory: '{run}/checkpoints'\n"
        "      checkpoint_capture_epochs:\n"+captures+
        f"      product_filename: '{run}/outputs/zhang_internal_products.csv'\n"
        f"      product_covariance_filename: '{run}/outputs/zhang_internal_product_covariance.csv'\n"
        f"outputs:\n  outputs_root: '{run}/outputs'\n"
        f"  metadata:\n    config_description: '{p['description']}'\n").encode()


def checkpoint_id(number,clock):
    return 'E29-'+clock.replace('-','').replace(' ','').replace(':','')+'-e'+str(number)+'-r32f626de9463'


def build_plan(raw,target_sha,auditor,profile):
    if profile not in PROFILES:raise ValueError('Only the two reviewed finite profiles are supported')
    if auditor!=AUDITOR:raise ValueError('Only the new isolated actual native auditor path')
    if len(target_sha)!=64 or any(c not in '0123456789abcdef' for c in target_sha) or target_sha==SOURCE_BINARY_SHA:
        raise ValueError('Distinct actual candidate SHA required')
    if len(raw)!=284226 or sha(raw)!=SOURCE_MANIFEST_SHA:raise ValueError('Pinned original source sidecar required')
    doc=json.loads(raw);source=P.decode(doc['config_text'].encode())
    if sha(doc['config_text'].encode())!=SOURCE_CONFIG_SHA or sha(doc['input_manifest_text'].encode())!=SOURCE_INPUT_SHA:
        raise ValueError('Original source canonical/input mismatch')
    if source['command_argument_count']!='149' or source['config_file_count']!='10' or doc['state_dimension']!=7274 or len(doc['sections'])!=9:
        raise ValueError('Original native inventory mismatch')
    f=PROFILES[profile];run=f['run'];overlay_path=run+'/config/'+f['overlay_name'];body=overlay(profile)
    target=dict(source);env={}
    for i in range(int(source['us20_global134_science_environment_count'])):
        prefix=f'us20_global134_science_environment_{i}_'
        env[target.pop(prefix+'key')]=target.pop(prefix+'value')
    env[FLAG]='1'
    env['ZHANG_R51_AR_START_GPST_SECONDS']=str(f['ar_start'])
    for i,(key,value) in enumerate(sorted(env.items())):
        prefix=f'us20_global134_science_environment_{i}_';target[prefix+'key']=key;target[prefix+'value']=value
    target['us20_global134_science_environment_count']=str(len(env))
    target.update(working_directory=run,us20_global134_producer_path=BINARY,
        us20_global134_checkpoint_policy=POLICY,effective_end_epoch=f['end'].replace(' ','T'),
        effective_max_epochs=str(f['max_epoch']),effective_outputs_root=run+'/outputs',
        effective_product_filename=run+'/outputs/zhang_internal_products.csv',
        effective_product_covariance_filename=run+'/outputs/zhang_internal_product_covariance.csv',
        effective_checkpoint_output_directory=run+'/checkpoints',config_file_count='11',command_argument_count='150')
    target['environment_ZHANG_R51_AR_START_GPST_SECONDS']=str(f['ar_start'])
    for suffix,value in dict(path=overlay_path,bytes=body.decode(),size=str(len(body)),sha256=sha(body)).items():
        target['config_file_10_'+suffix]=value
    original_argv=[source[f'command_argument_{i}'] for i in range(1,150)]
    argv=original_argv[:10]+[overlay_path]+original_argv[10:]
    for i,value in enumerate(argv,1):target[f'command_argument_{i}']=value
    for i in range(6):target.pop(f'effective_capture_epoch_{i}')
    target['effective_capture_epoch_count']=str(len(f['captures']))
    for i,(_,clock) in enumerate(f['captures']):target[f'effective_capture_epoch_{i}']=clock
    canonical=P.encode(P.actual_writer_order(source,target))
    env.update(ZHANG_US20_GLOBAL134_STAGE='AR',ZHANG_E29_RESTORE_BUNDLE=SOURCE,
        LD_LIBRARY_PATH=ROOT+'/candidate/lib')
    base=dict(argv=[BINARY]+argv,cwd=run,environment=env,binary_sha256=target_sha,launch_authorised_by_this_script=False)
    snapshot=run+'/restore_audit/ar_roundtrip_'+SOURCE_ID+'-validation-only'
    after=run+'/restore_audit/ar_roundtrip_'+SOURCE_ID+'-after-restore-validation-only'
    def restore(path):return dict(base,environment=dict(env,ZHANG_E29_RESTORE_BUNDLE=path,ZHANG_GEC134_RESTORE_ONLY='1'))
    def audit(before,after,report,kind,before_sha=target_sha):
        return dict(argv=[auditor,before,after,before_sha,target_sha,report,kind],required_argument_count=6)
    number,clock=f['captures'][-1];id=checkpoint_id(number,clock);own=run+'/checkpoints/'+id
    own_after=run+'/restore_audit/ar_roundtrip_'+id+'-validation-only'
    plan=dict(schema='CHECKPOINT_AUTHORITY_AR_PREPARE_V1',policy=POLICY,profile=profile,prepare_only=True,
        launch_approved=False,source_original_read_only=True,
        source=dict(directory=SOURCE,producer=SOURCE_BINARY,producer_sha256=SOURCE_BINARY_SHA,
            bundle_sha256=SOURCE_BUNDLE_SHA,manifest_sha256=SOURCE_MANIFEST_SHA,
            config_sha256=SOURCE_CONFIG_SHA,input_manifest_sha256=SOURCE_INPUT_SHA,
            completed_epoch=240,state_dimension=7274,sections=doc['sections']),
        target=dict(binary=BINARY,binary_sha256=target_sha,root=ROOT,run_root=run,
            new_epochs=f['max_epoch']-240,first_epoch=241,last_epoch=f['max_epoch'],
            start_gpst='2024-07-17 02:00:00',end_gpst=f['end'],
            captures=[dict(epoch=n,gpst=t) for n,t in f['captures']],AR_start=f['ar_start'],
            actual_AR_first_epoch=f['ar_first_epoch'],actual_AR_epochs=f['actual_ar_epochs'],
            canonical_sha256=sha(canonical),canonical_path='target_canonical_expected.txt'),
        overlay=dict(target_path=overlay_path,bytes=len(body),sha256=sha(body)),
        frozen_original_config_files=[dict(path=source[f'config_file_{i}_path'],bytes=int(source[f'config_file_{i}_size']),sha256=source[f'config_file_{i}_sha256']) for i in range(10)],
        ordered_pre_run_gates=[dict(kind='actual_PEA_source_restore_no_epoch',**restore(SOURCE)),
            dict(kind='independent_source_all9_33core',**audit(SOURCE,snapshot,run+'/restore_audit/source_e240_nine_section_independent.json','source-migration',SOURCE_BINARY_SHA)),
            dict(kind='actual_PEA_new_snapshot_restore_no_epoch',**restore(snapshot)),
            dict(kind='independent_snapshot_all9_33core',**audit(snapshot,after,run+'/restore_audit/new_e240_snapshot_nine_section_independent.json','same-producer'))],
        real_run_after_separate_root_gate=base,
        ordered_post_run_gates=[dict(kind='actual_PEA_own_checkpoint_restore_no_epoch',**restore(own)),
            dict(kind='independent_own_checkpoint_all9_33core',**audit(own,own_after,run+'/restore_audit/own_checkpoint_nine_section_independent.json','same-producer'))],
        acceptance=dict(required_complete_epochs=f['max_epoch']-240,required_complete_FLOAT_epochs=f['max_epoch']-240,
            required_actual_AR_epochs=f['actual_ar_epochs'],required_atomic_product_commits=f['max_epoch']-240,
            expected_float_epochs=f['max_epoch']-240,expected_AR_search_epochs=f['actual_ar_epochs'],
            expected_FLOAT_only_epochs=f['max_epoch']-240-f['actual_ar_epochs'],
            warmup_only_epochs=[],
            active_AR_epochs=[f['ar_first_epoch'],f['max_epoch']],AR_CALL_boundaries_are_not_scheduled_AR_count=True,
            native_exit0_is_not_sufficient=True,required_native_terminal_epoch=f['max_epoch'],
            source_format_export_raw_identity=True,own_AR_GEC_section_V2_must_be_proven=False,
            GEC_schema_must_match_actual_native_nonempty_registry=True,
            PPP_AR_source_schema_versions=[6,7],source_before_after_schema_must_match=True,
            receiver_exception_only='POST_EPOCH_READY_WORKSET_RESET_V1',all9_and33core_mandatory=True,
            ratio=2,joint_NIS=False,epoch_independent_AR=True,full_Q=True,WSL_GiB=20,
            raw_old_full_KF_e241_snapshot_available=False,scientific_accuracy_accepted=False),
        resource_policy=dict(new_resource_abort_thresholds=False,observational_samples_only=True))
    # Both half-hour checkpoints are actually restored in fresh processes,
    # with the same independent all-nine/raw33-core audit as the source gates.
    plan['ordered_post_run_gates']=[]
    for number,clock in f['captures']:
        identifier=checkpoint_id(number,clock);own=run+'/checkpoints/'+identifier
        own_after=run+'/restore_audit/ar_roundtrip_'+identifier+'-validation-only'
        plan['ordered_post_run_gates'] += [dict(kind='actual_PEA_own_checkpoint_restore_no_epoch',
            completed_epoch=number,**restore(own)),dict(kind='independent_own_checkpoint_all9_33core',
            **audit(own,own_after,run+'/restore_audit/own_e'+str(number)+'_nine_section_independent.json','same-producer'))]
    plan['acceptance']['old_algorithm_product_value_equality_required']=False
    plan['acceptance']['own_captured_checkpoint_actual_restores_required']=[n for n,_ in f['captures']]
    plan['acceptance']['formal_requires_passed_e300_replay']=True
    return plan,canonical,body


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--profile',choices=tuple(PROFILES),required=True)
    p.add_argument('--source-manifest',type=Path,required=True)
    p.add_argument('--target-binary-sha256',required=True)
    p.add_argument('--native-auditor',required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    args=p.parse_args()
    plan,canonical,body=build_plan(args.source_manifest.read_bytes(),args.target_binary_sha256,args.native_auditor,args.profile)
    if args.output_dir.exists():raise RuntimeError('Fresh versioned plan directory required')
    args.output_dir.mkdir(parents=True)
    for name,data in ((PROFILES[args.profile]['overlay_name'],body),('target_canonical_expected.txt',canonical),
        ('authority_hour_prepare_plan.json',(json.dumps(plan,indent=2)+'\n').encode())):
        with (args.output_dir/name).open('xb') as f:f.write(data)
    raw=(args.output_dir/'authority_hour_prepare_plan.json').read_bytes()
    print(json.dumps(dict(prepared=True,profile=args.profile,plan_sha256=sha(raw),
        no_compiler_or_solver_started=True),indent=2))


if __name__=='__main__':main()
