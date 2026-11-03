"""SHA-gated, single-use local AR hour and sequential real restore audits.

Preparing or importing never launches PEA. Each validation/formal/audit action
needs its own exact plan+nonce gate. No automatic retry, resource-abort monitor,
Windows publisher, server access, WSL restart or old-product value matching.
"""
from __future__ import annotations
import argparse
import ast
import csv
from datetime import datetime, timedelta, timezone
from functools import lru_cache
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import stat
import subprocess
import sys
import types

ART=Path('/mnt/c/Users/rx/Documents/GINAN/artifacts/checkpoint_authority_fix_20261104_01')
ROOT=Path('/home/rx/GINAN/checkpoint_authority_fix_20261104_01')
HERE=Path(__file__).resolve().parent
BASE_ART=ART.parent/'ar_p0_p4_5h_20261102_01'
BASE=BASE_ART/'harness/run_e240_5h_controller.py'
BASE_SHA='41c8aefee8453cce890c33af939b46481346c2f8b88d363be4f4d461a84a31ef'
CONTRACTS_SHA='03c67e08ebed4159c783fe169eccf60fb2e1d774d4ef19151b9d478c50173be3'
BASE_BUILDER_SHA='b78d170627eb44acc732fff9782113666cfa1dd6c63a19240a340fcec0a60125'
CODEC=ART.parent/'ar_checkpoint_codec_fix_20261103_01/harness/run_codec_replay_controller.py'
CODEC_SHA='553fbd479f2b2ca18d26e680425c372f5fb7864081ba9e8177b9b95b5c95c820'
OWNED=ART.parent/'us20_transport_goal_20261102_batch_01/e240_owned_runtime_v2.py'
OWNED_SHA='84bdff844634079f04463b7e1e2068d50331d925a79408609ee5b68f706406cc'
BUILDER=ART.parent/'datum_unknown_repair_20261103_01/harness/build_and_test_v4.py'
BUILDER_SHA='99cd029716a26c6bc186a253f3d1dd3b937aa954b833d0444874f8601cdc7a00'
SCHEMA='CHECKPOINT_AUTHORITY_AR_EXECUTION_PLAN_V1'
GATE_SCHEMA='CHECKPOINT_AUTHORITY_AR_ACTION_GATE_V1'
QUAL_SCHEMA='CHECKPOINT_AUTHORITY_CANDIDATE_QUALIFICATION_V1'
PREPARER_SHA='de31a7d0548cff80f47b4da2dcb01ad9a924ccf627080cd016a0a259c1ba3d3c'
ACTIONS=('validation','formal','audit')
PROFILES=('formal','replay_e300')
BASELINE_NATIVE=ART.parent/'datum_unknown_ar1h_20261103_01/build_01/NATIVE_RESULT.json'
BASELINE_NATIVE_SHA='f2d287e39899e10c034dc5764fddc215a1d8d9e08780b0d70de5444e1a67bfcb'
ADDITIONAL_NATIVE_CASES=('current_evidence_authority_snapshot','native_product_target_admission')
FAILED_RUN_ROOTS=(Path('/home/rx/GINAN/datum_unknown_ar1h_20261103_01/run'),
    Path('/home/rx/GINAN/datum_unknown_ar1h_20261103_01/controller_formal_01/formal'))
M=P=C=None
PROFILE=None


def require(ok,why):
    if not ok:raise RuntimeError(why)


def host_path(path):
    text=str(path).replace('\\','/')
    return Path(text[5]+':/'+text[7:]) if os.name=='nt' and text.startswith('/mnt/') else Path(text)


def sha(path):
    h=hashlib.sha256()
    with host_path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024**2),b''):h.update(block)
    return h.hexdigest()


def identity(path):
    path=host_path(path)
    return dict(path=str(path),bytes=path.stat().st_size,sha256=sha(path))


def strict_load(path):
    def pairs(items):
        result={}
        for k,v in items:
            require(k not in result,'Duplicate JSON key: '+k);result[k]=v
        return result
    def invalid(value):raise RuntimeError('Nonfinite JSON constant: '+value)
    def finite(value):
        number=float(value);require(math.isfinite(number),'Overflowed JSON number');return number
    return json.loads(host_path(path).read_text(encoding='utf-8'),object_pairs_hook=pairs,
        parse_constant=invalid,parse_float=finite)


def write_new(path,value):
    with Path(path).open('x',encoding='utf-8',newline='\n') as stream:
        json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n');stream.flush();os.fsync(stream.fileno())


def verify(row):
    require(isinstance(row,dict) and isinstance(row.get('path'),str) and type(row.get('bytes')) is int
        and row['bytes']>0 and re.fullmatch('[0-9a-f]{64}',row.get('sha256','')),'Malformed file pin')
    actual=identity(row['path'])
    require(actual['bytes']==row['bytes'] and actual['sha256']==row['sha256'],'Changed pinned file: '+row['path'])
    if 'resolved' in row:require(str(host_path(row['path']).resolve())==str(host_path(row['resolved']).resolve()),'Changed resolved file')


def import_pinned(path,digest,name):
    path=host_path(path);require(sha(path)==digest,'Changed reviewed helper: '+str(path))
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module;spec.loader.exec_module(module);return module


# Borrow only the reviewed pure parser/kernel-owner stage helpers. Rebinding
# their globals makes all roots, profiles and contracts below explicit. The
# old module's configure/prepare/launch functions are never invoked.
L=import_pinned(CODEC,CODEC_SHA,'_datum_frozen_codec_helpers')
for _name in ('proc_identity','terminal_stage','run_stage','marker_class','audit_epochs','emergency_failure'):
    _f=getattr(L,_name)
    globals()[_name]=types.FunctionType(_f.__code__,globals(),_name,_f.__defaults__,_f.__closure__)


@lru_cache(maxsize=2)
def epoch_grid(profile,resumed=False):
    require(profile in PROFILES and not resumed,'Only exact non-resumed finite profiles')
    count=120 if profile=='formal' else 60
    return tuple((datetime(2024,7,17,2)+timedelta(seconds=30*i)).strftime('%Y-%m-%d %H:%M:%S') for i in range(count))


def phase_for(profile,epoch):
    require(epoch in epoch_grid(profile),'Epoch outside exact finite grid');return 'AR'


def section_registry(manifest):
    expected={'GEC_AR_RUNTIME_V1':(1,2),'pea.controller.v1':(1,),
        'pea.receiver_runtime.v1':(3,),'pea.rinex_file_streams.v1':(2,),
        'pea.rng_declaration.v1':(1,),'pea.satellite_runtime.v1':(1,),
        'zhang_ambres_runtime':(3,),'zhang_graph_runtime':(5,),'zhang_ppp_ar_runtime':(6,7)}
    rows=manifest['sections'];indexed={row['name']:row for row in rows}
    require(len(rows)==len(indexed)==9 and set(indexed)==set(expected),'All exact nine sections required')
    for name,versions in expected.items():
        row=indexed[name]
        require(type(row['schema_version']) is int and row['schema_version'] in versions
            and type(row['payload_bytes']) is int and row['payload_bytes']>0
            and re.fullmatch('[0-9a-f]{64}',row['sha256']),'Invalid section framing: '+name)
    return indexed


def validate_nine_report(path,kind,before):
    report=strict_load(path);sections=section_registry(strict_load(Path(before)/'checkpoint_manifest.json'))
    require(report['passed'] is True and report['audit_kind']==kind and report['failed_fields']==[]
        and report['estimation_executed_by_auditor'] is False
        and report['whole_process_restore_must_be_proven_by_PEAs_native_receipt'] is True
        and report['scientific_accuracy_accepted'] is False,'Independent nine-section audit failed')
    checks=report['checks'];require(len([k for k in checks if k.startswith('core.')])==33
        and all(v is True for v in checks.values()),'All33 exact core checks required')
    indexed={row['name']:row for row in report['sections']}
    require(len(indexed)==len(report['sections'])==9 and set(indexed)==set(sections),'Nine independent sections absent')
    for name,row in indexed.items():
        require(type(row['schema']) is int and row['schema']==sections[name]['schema_version']
            and row['entire_persistent_payload_equal'] is True
            and row['declared_ready_workset_boundary_policy_only'] is (name=='pea.receiver_runtime.v1'),
            'Source-format whole section changed: '+name)
    require(report['state_dimension']>0 and report['P_entries_raw_IEEE']==report['state_dimension']**2,'Complete P absent')
    return identity(path)


def science_contract(science,directory):
    require(science.get('schema')=='CHECKPOINT_AUTHORITY_AR_PREPARE_V1' and science.get('profile')==PROFILE
        and science.get('prepare_only') is True and science.get('launch_approved') is False,'Wrong science authority')
    expected,canonical,overlay=P.build_plan(host_path(P.SOURCE+'/checkpoint_manifest.json').read_bytes(),
        science['target']['binary_sha256'],science['ordered_pre_run_gates'][1]['argv'][0],PROFILE)
    require(science==expected,'Science plan differs from exact finite renderer')
    directory=Path(directory);name=P.PROFILES[PROFILE]['overlay_name']
    require((directory/name).read_bytes()==overlay and
        (directory/'target_canonical_expected.txt').read_bytes()==canonical,'Canonical/overlay bytes changed')
    target=science['target'];wanted=epoch_grid(PROFILE)
    require(target['first_epoch']==241 and target['last_epoch']==240+len(wanted)
        and target['new_epochs']==len(wanted) and target['start_gpst']==wanted[0] and target['end_gpst']==wanted[-1],
        'Exact finite boundary changed')
    require(science['acceptance']['required_actual_AR_epochs']==len(wanted)
        and science['acceptance']['expected_FLOAT_only_epochs']==0
        and science['acceptance']['ratio']==2 and science['acceptance']['joint_NIS'] is False
        and science['acceptance']['epoch_independent_AR'] is True and science['acceptance']['WSL_GiB']==20
        and science['resource_policy']==dict(new_resource_abort_thresholds=False,observational_samples_only=True),
        'AR/resource policy changed')
    return directory/name,directory/'target_canonical_expected.txt'


def inventory(root):
    root=host_path(root)
    return {p.relative_to(root).as_posix():sha(p) for p in sorted((root/'src').rglob('*')) if p.is_file()}


def file_pins(node):
    if isinstance(node,dict):
        if {'path','bytes','sha256'}<=set(node):yield node
        else:
            for value in node.values():yield from file_pins(value)
    elif isinstance(node,list):
        for value in node:yield from file_pins(value)


def runtime_libraries(binary):
    result=subprocess.run(['/usr/bin/ldd',str(binary)],capture_output=True,text=True,check=True)
    require('not found' not in result.stdout,'Unresolved candidate runtime library')
    rows={}
    for line in result.stdout.splitlines():
        match=re.search(r'=>\s+(/\S+)',line) or re.match(r'\s*(/\S+)',line)
        if match:
            path=Path(match[1]).resolve(strict=True);rows[str(path)]=identity(path)
    require(rows,'No runtime libraries');return list(rows.values())


def failed_run_metadata(roots=FAILED_RUN_ROOTS):
    """No input traversal or large-log rehash: preserve exact lstat identities."""
    result={}
    for root in roots:
        root=Path(root)
        require(root.is_dir() and not root.is_symlink() and root.resolve(strict=True)==root,
            'Original failed run root disappeared or was redirected')
        entries={}
        for directory,subdirs,files in os.walk(root,followlinks=False):
            for path in [Path(directory),*[Path(directory)/name for name in sorted(subdirs+files)]]:
                key=path.relative_to(root).as_posix()
                if key in entries:continue
                info=path.lstat()
                kind='link' if stat.S_ISLNK(info.st_mode) else 'directory' if stat.S_ISDIR(info.st_mode) else 'file'
                require(kind!='file' or stat.S_ISREG(info.st_mode),'Unexpected failed-run special file')
                row=dict(kind=kind,device=info.st_dev,inode=info.st_ino,bytes=info.st_size,mtime_ns=info.st_mtime_ns)
                if kind=='link':row['target']=os.readlink(path)
                entries[key]=row
        result[str(root)]=entries
    return result


def qualification(path,digest):
    require(sha(path)==digest,'Exact candidate qualification SHA required');q=strict_load(path)
    require(q.get('schema')==QUAL_SCHEMA and q.get('passed') is True,'Wrong qualification')
    for k in ('candidate','auditor','build_plan','build_result','native_result'):verify(q[k])
    require(q['candidate']['path']==str(ROOT/'candidate/pea') and
        q['auditor']['path']==str(ROOT/'candidate/us20_e240_nine_section_audit'),'Unexpected deployed executable path')
    bp,br,nr=(strict_load(q[k]['path']) for k in ('build_plan','build_result','native_result'))
    require(bp['schema']=='DATUM_UNKNOWN_ISOLATED_BUILD_V4' and bp['ratio']==2 and bp['joint_nis'] is False
        and bp['epoch_independent_ar'] is True and bp['solver_run'] is False and br['passed'] is True
        and br['plan_sha256']==q['build_plan']['sha256'] and nr['plan_sha256']==q['build_plan']['sha256'],
        'Build/native plan identity differs')
    require(len(bp['compile'])==len(br['compile'])==158 and len({r['object'] for r in br['compile']})==158
        and {(r['source'],r['object']) for r in bp['compile']}=={(r['source'],r['object']) for r in br['compile']}
        and all(r['passed'] is True and sha(r['object'])==r['object_sha256'] for r in br['compile']),
        'Full158 production object evidence missing')
    require(br['link']['passed'] is True and br['link']['exit_code']==0
        and br['link']['command']==bp['full_link_command']
        and {x for x in br['link']['command'] if x.endswith('.o')}=={r['object'] for r in br['compile']}
        and sha(br['candidate'])==br['candidate_sha256']==q['candidate']['sha256'],
        'Actual full158 candidate link differs')
    require(inventory(bp['source'])==bp['source_manifest']==inventory(q['source_input']),
        'Full frozen source/header tree differs')
    for row in bp['compile']:
        require(isinstance(row['project_dependencies'],list) and
            all(isinstance(k,str) and k in bp['source_manifest'] for k in row['project_dependencies']),
            'Compile dependencies not bound to frozen source')
    require(sha(BASELINE_NATIVE)==BASELINE_NATIVE_SHA,'Baseline complete regression inventory changed')
    baseline_native=strict_load(BASELINE_NATIVE)
    require(baseline_native['passed'] is True and len(baseline_native['tasks'])==51,
        'Baseline regression inventory is not the complete actual suite')
    wanted_names={r['name'] for r in baseline_native['tasks']}|set(ADDITIONAL_NATIVE_CASES)
    planned_names={r['name'] for r in bp['native_tasks']}
    actual_names={r['name'] for r in nr['tasks']}
    require(len(bp['native_tasks'])==len(planned_names)==len(nr['tasks'])==len(actual_names)==len(wanted_names)
        and planned_names==actual_names==wanted_names,'Complete baseline plus new regression suite required')
    execute_count=sum(r['execute'] is True for r in bp['native_tasks'])
    require(execute_count==len(wanted_names)-1 and all(r['execute'] is True or
        r['name']=='us20_e240_nine_section_audit' for r in bp['native_tasks']),
        'Only the independent auditor may be build-only')
    require(nr['schema']=='DATUM_UNKNOWN_NATIVE_AGGREGATE_V4' and nr['passed'] is True
        and nr['production_sha256']==q['candidate']['sha256'] and nr['run_cases']==execute_count and nr['build_only_cases']==1
        and nr['full_production_link'] is True and nr['actual_executions_in_this_attempt']==execute_count
        and nr['accepted_earlier_current_candidate_executions']==0 and nr['accepted_earlier_build_only_records']==0
        and all(r['passed'] is True for r in nr['tasks']),
        'All current-producer native contracts must actually execute in this attempt')
    for name in ('gec_ar_runtime_checkpoint','gps_prepared_support_contract','historical_evidence_archive',
                 'us20_e240_continuation_policy',*ADDITIONAL_NATIVE_CASES):
        rows=[r for r in nr['tasks'] if r['name']==name]
        require(len(rows)==1 and rows[0]['execute_requested'] is True and rows[0]['execute']['passed'] is True
            and rows[0]['execute']['exit_code']==0,'Mandatory actual native execution missing: '+name)
    auditors=[r for r in nr['tasks'] if r['name']=='us20_e240_nine_section_audit']
    require(len(auditors)==1 and auditors[0]['execute_requested'] is False and
        auditors[0]['executable']['sha256']==q['auditor']['sha256'],
        'Independent auditor is not the actual built driver')
    snapshot_case=next(r for r in nr['tasks'] if r['name']=='current_evidence_authority_snapshot')
    require('CURRENT_EVIDENCE_AUTHORITY_SNAPSHOT_PASS' in host_path(snapshot_case['execute']['log']).read_text(),
        'New authority rollback regression did not report its actual execution')
    admission_case=next(r for r in nr['tasks'] if r['name']=='native_product_target_admission')
    require('NATIVE_PRODUCT_TARGET_ADMISSION_PASS' in host_path(admission_case['execute']['log']).read_text(),
        'New pending-target regression did not report its actual execution')
    pins=[identity(BASELINE_NATIVE)]+list(file_pins(q))+list(file_pins(nr))+list(file_pins(bp.get('tool_and_link_pins',{})))
    for path,value in bp['external_dependencies'].items():
        pin=identity(path);require(pin['sha256']==value,'External build dependency changed: '+path);pins.append(pin)
    return q,bp,br,nr,pins


def prepare_execution(args):
    M.local_only();require(not M.EVIDENCE.exists() and not M.CONTROL.exists() and not M.RUN.exists(),
        'Fresh canonical plan/run/controller required; no automatic retry')
    require(ROOT.is_dir() and not ROOT.is_symlink() and ROOT.resolve(strict=True)==ROOT,'Canonical isolated root required')
    require(sha(args.pure_test_receipt)==args.pure_test_sha256,'Exact pure-test receipt SHA required')
    pure=strict_load(args.pure_test_receipt)
    require(pure.get('passed') is True and pure.get('solver_processes_started')==0 and pure.get('tests_run',0)>=29
        and pure['controller']['sha256']==sha(__file__) and pure['test_source']['sha256']==sha(HERE/'test_authority_hour.py')
        and pure['preparer']['sha256']==sha(HERE/'prepare_authority_hour.py'),'Current pure source/test receipt required')
    q,bp,br,nr,pins=qualification(args.qualification,args.qualification_sha256)
    M.pin_argument(args.science_plan,args.science_sha256);science=strict_load(args.science_plan)
    overlay,canonical=science_contract(science,Path(args.science_plan).parent)
    require(science['target']['binary_sha256']==q['candidate']['sha256'] and
        science['ordered_pre_run_gates'][1]['argv'][0]==q['auditor']['path'],'Science executables not qualified')
    source=Path(P.SOURCE);doc=strict_load(source/'checkpoint_manifest.json')
    require(sha(source/'checkpoint_manifest.json')==P.SOURCE_MANIFEST_SHA and
        sha(source/'checkpoint.bundle')==P.SOURCE_BUNDLE_SHA,'Original e240 exact bundle changed')
    inputs=C.input_records(doc['input_manifest_text']);C.verify_inputs(inputs)
    configs=science['frozen_original_config_files'];require(len(configs)==10,'All10 original configurations required')
    for row in configs:verify(row)
    helpers=M.helper_pins();prerequisites=strict_load(M.HELPERS/'native_prerequisites.json')
    reader,receiver=prerequisites['independent_core_reader'],prerequisites['independent_receiver_auditor']
    verify(reader);verify(receiver)
    old_art=ART.parent/'us20_global134_osb5h_20261031_01'
    old_bp,old_br=(strict_load(old_art/name) for name in ('BUILD_PLAN.json','BUILD_RESULT.json'))
    original='/home/rx/GINAN/us20_global134_osb5h_20261031_v1/source'
    require(old_br['passed'] is True and old_br['candidate']['sha256']==P.SOURCE_BINARY_SHA
        and inventory(original)==old_bp['source_manifest'],'Original producer source changed')
    original_objects=[identity(r['object']) for r in old_br['compile']]
    require(len(original_objects)==158 and all(x['sha256']==r['object_sha256'] for x,r in zip(original_objects,old_br['compile'])),
        'Original158 object tree changed')
    preserved=M.failed_validation_pins()
    pins += helpers+configs+original_objects+[identity(p) for p in
        (__file__,HERE/'prepare_authority_hour.py',HERE/'test_authority_hour.py',args.pure_test_receipt,args.qualification,
         BASE,BASE.parent/'e240_5h_controller_contracts.py',BASE.parent/'build_5h_candidate.py',
         BASE.parent/'prepare_e240_5h_plan.py',CODEC,OWNED,BUILDER,source/'checkpoint.bundle',source/'checkpoint_manifest.json',
         P.SOURCE_BINARY,old_art/'BUILD_PLAN.json',old_art/'BUILD_RESULT.json',overlay,canonical,'/usr/bin/python3','/usr/bin/time')]
    unique={}
    for row in pins:
        pin={k:row[k] for k in ('path','bytes','sha256')}
        require(pin['path'] not in unique or unique[pin['path']]==pin,'Conflicting file identity')
        unique[pin['path']]=pin
    for pin in unique.values():verify(pin)
    runtime=runtime_libraries(ROOT/'candidate/pea')
    M.EVIDENCE.mkdir();M.CONTROL.mkdir()
    plan=dict(schema=SCHEMA,science_plan=identity(args.science_plan),science=science,
        qualification=identity(args.qualification),final_source_manifest=bp['source_manifest'],source_copy=bp['source'],
        source_input=q['source_input'],original_source=original,original_source_manifest=old_bp['source_manifest'],
        original_objects=original_objects,input_records=inputs,pins=list(unique.values()),
        failed_run_lstat_identity=failed_run_metadata(),
        executables=[q['candidate'],q['auditor'],reader,receiver],runtime_libraries=runtime,
        runtime_lib_directory=bp['runtime_lib_directory'],helper_pins=helpers,canonical=identity(canonical),overlay=identity(overlay),
        cwd_relative_input_bindings=[M.relative_input_binding(science)],preserved_failed_validation_evidence=preserved,
        action_gates_are_separate=True,no_solver_started_by_prepare=True,no_resource_abort_thresholds=True,
        no_windows_publisher=True,no_automatic_retry=True,scientific_accuracy_accepted=False)
    write_new(M.EVIDENCE/'EXECUTION_PLAN.json',plan)
    print(json.dumps(dict(prepared=True,plan=identity(M.EVIDENCE/'EXECUTION_PLAN.json'),solver_started=False),indent=2))


def verify_plan(path,digest,full_inputs=False):
    M.pin_argument(path,digest);plan=strict_load(path)
    require(plan['schema']==SCHEMA and plan['no_resource_abort_thresholds'] is True
        and plan['no_windows_publisher'] is True and plan['no_automatic_retry'] is True,'Wrong execution policy')
    require(Path(path).resolve()==(M.EVIDENCE/'EXECUTION_PLAN.json').resolve(),'Unexpected plan path')
    verify(plan['science_plan']);science_contract(plan['science'],Path(plan['science_plan']['path']).parent)
    M.verify_relative_input_bindings(plan)
    require(plan['preserved_failed_validation_evidence']==M.failed_validation_pins(),'Old failed evidence changed')
    require(plan['failed_run_lstat_identity']==failed_run_metadata(),
        'Original failed run files/links or metadata changed')
    for row in plan['pins']+plan['executables']:verify(row)
    require(inventory(plan['source_copy'])==plan['final_source_manifest']==inventory(plan['source_input'])
        and inventory(plan['original_source'])==plan['original_source_manifest'],'Current/original source tree changed')
    for row in plan['original_objects']:verify(row)
    require(runtime_libraries(ROOT/'candidate/pea')==plan['runtime_libraries'],'Candidate dynamic library set changed')
    require((ROOT/'candidate/lib').is_symlink() and (ROOT/'candidate/lib').resolve()==Path(plan['runtime_lib_directory']).resolve(),
        'Candidate library directory binding changed')
    if full_inputs:C.verify_inputs(plan['input_records'])
    return plan


def replay_acceptance(pin):
    expected=ART/'controller_replay_e300_01/audit_RESULT.json'
    require(isinstance(pin,dict) and pin.get('path')==str(expected),
        'Formal execution requires this candidate\'s actual e300 replay audit')
    verify(pin);result=strict_load(expected)
    replay_plan=ART/'controller_replay_e300_01/EXECUTION_PLAN.json'
    formal_plan=ART/'controller_formal_01/EXECUTION_PLAN.json'
    replay,formal=strict_load(replay_plan),strict_load(formal_plan)
    require(result.get('passed') is True and result.get('profile')=='replay_e300'
        and result.get('execution_plan_sha256')==sha(replay_plan)
        and result.get('complete_epochs')==result.get('atomic_product_commits')==60
        and result.get('complete_FLOAT_epochs')==result.get('AR_enabled_epochs')==60
        and result.get('all_captured_checkpoints_actual_roundtrip_tested') is True
        and len(result.get('actual_checkpoint_restores',[]))==1
        and len(result.get('products',[]))==60
        and len(result.get('regression_product_events',[]))==5
        and set(result.get('checkpoints',{}))=={'300'}
        and result.get('original_source_objects_configuration_inputs_unchanged') is True,
        'Actual e300 replay, fresh restore and independent full audit must pass before formal execution')
    require(replay['qualification']==formal['qualification']
        and replay['science']['source']==formal['science']['source']
        and replay['science']['target']['binary_sha256']==formal['science']['target']['binary_sha256'],
        'Replay proof belongs to a different candidate or original source')
    return pin


def validate_gate(action,gate,execution_sha):
    require(gate.get('schema')==GATE_SCHEMA and gate.get('action')==action and gate.get('profile')==PROFILE
        and gate.get('approved') is True and gate.get('execution_plan_sha256')==execution_sha
        and gate.get('automatic_retry_authorised') is False,'Exact separate root action gate required')
    if PROFILE=='formal':replay_acceptance(gate.get('replay_acceptance'))
    if action in ('formal','audit'):
        previous=M.action_files('validation' if action=='formal' else 'formal')['result']
        require(gate.get('predecessor_result')==identity(previous) and strict_load(previous)['passed'] is True
            and strict_load(previous)['execution_plan_sha256']==execution_sha,'Passed exact predecessor required')
        result=strict_load(previous)
        if action=='formal':
            require(result['actual_fresh_PEA_restores']==2 and result['all_nine_section_audits']==2
                and len(result['stages'])==2 and result['finite_data_run_started'] is False,
                'Both real source restores and independent audits must precede formal AR')
        else:
            require(result['complete_epochs']==result['complete_FLOAT_epochs']==result['AR_enabled_epochs']==len(epoch_grid(PROFILE)),
                'All complete actual epochs must precede final audit')


def restore_and_audit(plan,owned,label,stage,audit,completed):
    native,_,markers=run_stage(owned,label,stage)
    require(markers.restore_seen and markers.strict_restore_seen and not markers.epochs,'Fresh restore consumed an epoch')
    before,after=Path(audit['argv'][1]),Path(audit['argv'][2])
    receipt=M.native_restore_receipt(plan,before,after,completed)
    source_sections=section_registry(strict_load(before/'checkpoint_manifest.json'))
    after_sections=section_registry(strict_load(after/'checkpoint_manifest.json'))
    require(all(source_sections[n]['schema_version']==after_sections[n]['schema_version'] for n in source_sections),
        'Source-format restore changed a section version')
    receiver=M.independent_receiver(plan,owned,after,label,stage['environment']);terminal_stage(receiver['independent_process'])
    independent,_,_=run_stage(owned,label+'_nine_sections',dict(argv=audit['argv'],cwd=stage['cwd'],environment=stage['environment']))
    report=validate_nine_report(audit['argv'][5],audit['argv'][6],before)
    return dict(restore=native,native=receipt,receiver=receiver,independent_audit=independent,independent_report=report)


def run_validation(plan,owned):
    if PROFILE=='replay_e300':
        M.RUN.parent.mkdir(exist_ok=True)
        require(M.RUN.parent.is_dir() and not M.RUN.parent.is_symlink() and M.RUN.parent.resolve(strict=True)==M.RUN.parent,
            'Isolated validation parent redirected')
    M.materialize(plan);M.verify_relative_input_bindings(plan,require_materialized=True)
    stages=plan['science']['ordered_pre_run_gates'];results=[]
    for i in (0,2):results.append(restore_and_audit(plan,owned,'source_restore_'+str(i//2+1),stages[i],stages[i+1],240))
    return dict(passed=True,profile=PROFILE,actual_fresh_PEA_restores=2,all_nine_section_audits=2,
        stages=results,finite_data_run_started=False)


def run_formal(plan,owned):
    M.verify_relative_input_bindings(plan,require_materialized=True)
    require(not list((M.RUN/'outputs').glob('*.epochs/*/COMMIT.json')),'Empty validated product sink required')
    native,out,markers=run_stage(owned,'finite_'+PROFILE,plan['science']['real_run_after_separate_root_gate'])
    rows=audit_epochs(markers,PROFILE);terminal=240+len(epoch_grid(PROFILE))
    tail,evidence=C.bounded_tail(out/'runner.log')
    require(f'Exiting at epoch {terminal} as epoch count {terminal} has been reached' in tail,'Natural exact finite terminal missing')
    return dict(passed=True,profile=PROFILE,native=native,epochs=rows,complete_epochs=len(rows),complete_FLOAT_epochs=len(rows),
        AR_CALL_timing_envelopes=len(rows),AR_enabled_epochs=len(rows),
        epochs_with_R51_SEARCH_ENTERED=sum(bool(r.get('actual_search_marker_seen')) for r in rows),
        AR_epoch_count_semantics='AR_ENABLED transactions, not guaranteed LAMBDA calls',
        terminal_barrier_bounded_tail=evidence,all_products_and_checkpoints_independently_audited=False)


def validate_event_product(commit,rows):
    require(commit.get('route')=='ENHANCED' and commit.get('ar') is True,
        'Original failure event regressed to a global FLOAT_ONLY publication')
    systems=commit.get('systems',{})
    require(set(systems)=={'GPS','GAL','BDS'},'Three-system event publication is incomplete')
    usable={}
    for system in ('GPS','GAL','BDS'):
        count=systems[system].get('pppar_usable_phase_rows')
        require(type(count) is int and count>=0 and count%2==0,'Unpaired usable phase count')
        usable[system]=count//2
    require(usable['GAL']>=2 and usable['BDS']>=2,
        'Galileo/BDS lost their valid two-frequency products during the GPS pending-target event')
    require(usable['GPS']==0 or usable['GPS']>=2,'A single GPS target must remain pending')
    pending_rows=[]
    if usable['GPS']==0:
        pending_rows=[row for row in rows if row.get('satellite','').startswith('G')]
        for row in pending_rows:
            require(all(row.get(field)=='0' for field in ('ar_valid','pppar_usable','integer_valid',
                'integer_transport_valid','dual_frequency_ar_valid')),
                'Pending GPS target retained an AR/integer publication grant')
    return dict(epoch=commit['epoch'],passed=True,route=commit['route'],dual_usable=usable,
        pending_GPS_rows=len(pending_rows),pending_GPS_rows_have_no_AR_or_integer_grant=True)


def regression_product_events(paths):
    expected=[(datetime(2024,7,17,2,27,30)+timedelta(seconds=30*i)).strftime('%Y-%m-%d %H:%M:%S') for i in range(5)]
    indexed={}
    for path in paths:
        commit=strict_load(path)
        if commit.get('epoch') in expected:
            require(commit['epoch'] not in indexed,'Duplicate regression event epoch')
            with (path.parent/'products.csv').open(newline='') as stream:rows=list(csv.DictReader(stream))
            indexed[commit['epoch']]=validate_event_product(commit,rows)
    require(sorted(indexed)==expected,'All five original failure epochs must be replayed and independently audited')
    return [indexed[epoch] for epoch in expected]


def run_audit(plan,owned):
    M.verify_relative_input_bindings(plan,require_materialized=True)
    auditor=M.imported(M.HELPERS/'product_auditor.py','datum_atomic_product_auditor')
    paths=sorted((M.RUN/'outputs/zhang_internal_products.csv.epochs').glob('*/COMMIT.json'))
    reports=[auditor.audit_bundle(p.parent) for p in paths]
    require(len(reports)==len(epoch_grid(PROFILE)) and sorted(r['epoch'] for r in reports)==list(epoch_grid(PROFILE)),
        'Missing/extra/duplicate atomic full-Q product bundle')
    regression=regression_product_events(paths)
    commits=[strict_load(path) for path in paths]
    commit_by_epoch={row['epoch']:row for row in commits}
    captures={};restores=[];gates=plan['science']['ordered_post_run_gates']
    for index,capture in enumerate(plan['science']['target']['captures']):
        number=capture['epoch'];identifier=P.checkpoint_id(number,capture['gpst']);checkpoint=M.RUN/'checkpoints'/identifier
        captures[str(number)]=M.content_audit(plan,owned,checkpoint,number)
        sections=section_registry(strict_load(checkpoint/'checkpoint_manifest.json'))
        require(sections['zhang_ppp_ar_runtime']['schema_version']==7,'New actual AR checkpoint must retain PPP predictor schema7')
        publication=commit_by_epoch[capture['gpst']]
        active_products=publication.get('ar') is True and any(
            value.get('pppar_usable_phase_rows',0)>0 for value in publication.get('systems',{}).values())
        if active_products:
            require(sections['GEC_AR_RUNTIME_V1']['schema_version']==2,
                'Active multi-system products cannot be checkpointed as an empty legacy GEC section')
        # GEC V1 is permitted only when it really is empty; the native exporter
        # and full independent raw audit remain the authority for its contents.
        restores.append(restore_and_audit(plan,owned,'own_e'+str(number),gates[2*index],gates[2*index+1],number))
    return dict(passed=True,profile=PROFILE,complete_epochs=len(reports),atomic_product_commits=len(reports),
        complete_FLOAT_epochs=len(reports),AR_enabled_epochs=len(reports),
        epochs_with_R51_SEARCH_ENTERED=strict_load(M.action_files('formal')['result'])['epochs_with_R51_SEARCH_ENTERED'],
        products=reports,checkpoints=captures,actual_checkpoint_restores=restores,
        regression_product_events=regression,
        enhanced_AR_product_epochs=sum(row.get('route')=='ENHANCED' and row.get('ar') is True for row in commits),
        float_only_product_epochs=sum(row.get('route')=='FLOAT_ONLY' for row in commits),
        all_captured_checkpoints_actual_roundtrip_tested=True,
        scientific_accuracy_accepted=False,old_algorithm_product_equality_performed=False)


def configure(profile):
    global M,P,C,PROFILE
    require(profile in PROFILES,'Only the two exact finite profiles');PROFILE=profile
    Ppath=host_path(HERE/'prepare_authority_hour.py')
    require(sha(Ppath)==PREPARER_SHA,'Changed exact finite preparer')
    spec=importlib.util.spec_from_file_location('_datum_hour_preparer',Ppath);P=importlib.util.module_from_spec(spec);spec.loader.exec_module(P)
    for path,digest in ((BASE,BASE_SHA),(BASE.parent/'e240_5h_controller_contracts.py',CONTRACTS_SHA),
        (BASE.parent/'build_5h_candidate.py',BASE_BUILDER_SHA),(BASE.parent/'prepare_e240_5h_plan.py',P.BASE_SHA),
        (OWNED,OWNED_SHA),(BUILDER,BUILDER_SHA)):
        require(sha(path)==digest,'Changed lifecycle dependency: '+str(path))
    tree=ast.parse(host_path(BASE).read_text(encoding='utf-8'))
    mapping={'gates_in_order':'ordered_pre_run_gates',
        'formal_run_after_root_review_and_all_gates_pass':'real_run_after_separate_root_gate',
        'gates_after_formal_run_without_changing_original_product_sinks':'ordered_post_run_gates',
        'five_hour_result.json':'hour_result.json','five_hour_complete':'hour_complete',
        'US20_E240_ACTION_RESULT_V1':'CHECKPOINT_AUTHORITY_AR_ACTION_RESULT_V1'}
    class Adapt(ast.NodeTransformer):
        def visit_Constant(self,node):
            if isinstance(node.value,str) and node.value in mapping:return ast.copy_location(ast.Constant(mapping[node.value]),node)
            return node
        def visit_FunctionDef(self,node):
            node=self.generic_visit(node)
            if node.name=='run':
                class Counts(ast.NodeTransformer):
                    def visit_Constant(self,n):
                        if type(n.value) is int and n.value==600:return ast.copy_location(ast.parse('len(C.epochs())',mode='eval').body,n)
                        if type(n.value) is int and n.value==840:return ast.copy_location(ast.parse('240+len(C.epochs())',mode='eval').body,n)
                        return n
                    def visit_keyword(self,n):
                        n=self.generic_visit(n)
                        if n.arg=='five_hour_complete':n.arg='hour_complete'
                        if n.arg=='actual_fresh_PEA_restores_after_formal':n.value=ast.parse('len(C.capture_schedule())',mode='eval').body
                        if n.arg=='checkpoints_not_actual_roundtrip_tested':n.value=ast.List(elts=[],ctx=ast.Load())
                        if n.arg=='all_ten_checkpoints_actual_roundtrip_tested':
                            n.arg='all_captured_checkpoints_actual_roundtrip_tested';n.value=ast.Constant(True)
                        return n
                return Counts().visit(node)
            return node
    module=types.ModuleType('_datum_owned_lifecycle');module.__file__=str(host_path(BASE));sys.modules[module.__name__]=module
    exec(compile(ast.fix_missing_locations(Adapt().visit(tree)),str(BASE),'exec'),module.__dict__);M=module;C=M.C
    C.manifest_sections=section_registry;C.Markers=marker_class(C,profile);C.epochs=lambda:list(epoch_grid(profile))
    C.capture_schedule=lambda:[dict(number=n,gpst=t,clock=t[11:].replace(':','')) for n,t in P.PROFILES[profile]['captures']]
    def checkpoint_epoch(number):
        rows=[r['gpst'] for r in C.capture_schedule() if r['number']==number];require(len(rows)==1,'Unknown capture');return rows[0]
    C.checkpoint_epoch=checkpoint_epoch
    p=types.SimpleNamespace(**vars(P));p.RUN=P.PROFILES[profile]['run'];p.SNAPSHOT=p.RUN+'/restore_audit/ar_roundtrip_'+P.SOURCE_ID+'-validation-only'
    p.OVERLAY_PATH=p.RUN+'/config/'+P.PROFILES[profile]['overlay_name']
    M.ART,M.ROOT,M.RUN,M.P=ART,ROOT,Path(p.RUN),p
    M.CONTROL=ROOT/('controller_'+profile+'_01');M.EVIDENCE=ART/('controller_'+profile+'_01')
    M.SCRIPT,M.SCHEMA,M.OWNED_RUNTIME=Path(__file__).resolve(),SCHEMA,OWNED
    M.prepare,M.verify_plan,M.validate_gate=prepare_execution,verify_plan,validate_gate
    M.science_contract=science_contract;M.run_validation,M.run_formal,M.run_audit=run_validation,run_formal,run_audit
    M.load,M.write_new,M.emergency_startup_failure=strict_load,write_new,emergency_failure
    return M


def main():
    parser=argparse.ArgumentParser(description=__doc__);mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare',action='store_true');mode.add_argument('--launch',choices=ACTIONS);mode.add_argument('--run',choices=ACTIONS)
    parser.add_argument('--profile',choices=PROFILES)
    for name in ('science-plan','execution-plan','gate','qualification','pure-test-receipt'):parser.add_argument('--'+name,type=Path)
    for name in ('science-sha256','plan-sha256','gate-sha256','qualification-sha256','pure-test-sha256','action-nonce'):
        parser.add_argument('--'+name)
    args=parser.parse_args();profile=args.profile
    if profile is None:
        require(not args.prepare and args.execution_plan and sha(args.execution_plan)==args.plan_sha256,'Child exact plan required')
        profile=strict_load(args.execution_plan)['science']['profile']
    engine=configure(profile)
    if args.prepare:engine.prepare(args)
    elif args.launch:engine.launch(args)
    else:engine.run(args)


if __name__=='__main__':main()
