"""Isolated serial authority-repair production and full native build.

Reuse only SHA-verified actual dependency-closed production objects from the
8900 build. Every native case is freshly compiled, linked and executed (except
the independent audit driver, which is built here and actually run on archives
by the separate owned workflow). Old failed runs and evidence remain unchanged.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import time
import re

ART = Path('/mnt/c/Users/rx/Documents/GINAN/artifacts/checkpoint_authority_fix_20261104_01')
REPAIR_ART = ART.parent/'datum_unknown_repair_20261103_01'
ROOT = Path('/home/rx/GINAN/checkpoint_authority_fix_20261104_01')
INPUT = ART/'worktree/source_snapshots/global134_and_user_20261103/server'
BASE = ART.parent/'datum_unknown_ar1h_20261103_01/build_01'
PLAN_SHA = '86f5a49584336fec9bf94c36c98aee4e423b713071470c8a5505bb4f8cf3642e'
RESULT_SHA = 'a6184db76622f067d4b6b21800d4bd2a87e2f57e913c769dfe7e8635a2a64c26'
NATIVE_BASE = ART.parent/'ar_checkpoint_codec_fix_20261103_01/build_attempt_03'
NATIVE_PLAN_SHA = '9d821b96c7bc4b303204ecf96b291b4624293024f8f841e940e3d4847925ad71'
NATIVE_RESULT_SHA = '5db2d148e383679a2580e0b23dd844a3d91e41fabd1b36c5db13c502cbe4a9bc'
NATIVE_RECIPE_SHA = '52d0e26ec94f62e1b7e959e6628fd621e6179f1385323dbe7c06aa7f4577a4f4'
PRIOR_NATIVE = REPAIR_ART/'build_04'
PRIOR_PLAN_SHA = '7cb24438a374e3064d95c13b43a56a4a3e2b0d870453e295c4d9bf9fb0fdfde2'
PRIOR_BUILD_SHA = 'd42254f1d2fe810bedd414fdc8b8c60a01819532a37cf1c87e72c78bcbf58abf'
PRIOR_NATIVE_SHA = 'aa1d07b4eb7c23433a09db0e1e83c86be5d8bd43f5ba310600829fbe5bb62f28'
ANCESTOR_NATIVE = REPAIR_ART/'build_03'
ANCESTOR_RECEIPTS = {'BUILD_PLAN.json':'3350dd1e7c85f5877148ff6ae5885ce2e266231950f097c2abf2766cb642cd51',
    'BUILD_RESULT.json':'792f11383ab19070035d1b0afb5e0c749dcc4836e72f2d08e536c1216a659440',
    'NATIVE_RESULT.json':'432468f69ab828a1e8eae7411d23ac68411e05005381fa1427857eececdf68b1'}
PRIOR_PRODUCER_SHA = '00a16f65eee965038b672b5ad773468b5798d86b1690cecf1771902b6474f657'
SUPPORT_FILES = {'gec_private_ambres_wire_schema.hpp','generate_codec_test_schema.py','PROVENANCE.json'}
def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1048576),b''): h.update(b)
    return h.hexdigest()
def load(path): return json.loads(Path(path).read_text())
def save(path, data):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_name(path.name+'.tmp')
    temp.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n'); os.replace(temp,path)
def require(ok, why):
    if not ok: raise RuntimeError(why)
def inventory(root):
    return {p.relative_to(root).as_posix():sha(p) for p in sorted((root/'src').rglob('*')) if p.is_file()}
def support_inventory(root):
    return {p.relative_to(root).as_posix():sha(p) for p in sorted((root/'test_support').rglob('*')) if p.is_file()}
def file_pin(path):
    path=Path(path)
    return dict(path=str(path),resolved=str(path.resolve()),bytes=path.stat().st_size,sha256=sha(path))
def verify_pin(pin):
    require(Path(pin['path']).is_file() and sha(pin['path'])==pin['sha256'],'Pinned file changed '+pin['path'])
    if 'resolved' in pin: require(str(Path(pin['path']).resolve())==pin['resolved'],'Pinned symlink changed '+pin['path'])
def verify_support(root, manifest):
    require(support_inventory(root)==manifest,'Frozen test_support changed')
    generator=root/'test_support/generate_codec_test_schema.py'
    namespace={'__name__':'frozen_codec_schema_generator'}
    exec(compile(generator.read_text(),str(generator),'exec'),namespace)
    companion=root/'src/tests/run_gec_ar_runtime_checkpoint.py'
    header,proof=namespace['render'](root,sha(companion))
    require(header.encode()==(root/'test_support/gec_private_ambres_wire_schema.hpp').read_bytes(),
            'Private codec declaration mirror is not generated from current frozen producer')
    original=load(root/'test_support/PROVENANCE.json')
    for name in ('source_sha256','companion_sha256','extraction_function_sha256','declaration_names','header_sha256'):
        require(original[name]==proof[name],'Codec declaration provenance differs from frozen producer '+name)
    require(original['validators_or_exporters_copied'] is False,'Codec mirror must contain declarations only')
    return proof
def dep_tokens(path):
    text=Path(path).read_text().replace('\\\n',' ')
    return [os.path.normpath(p) for p in shlex.split(text.split(':',1)[1])]
def source_root(command):
    src=command[command.index('-c')+1]
    return src.split('/src/',1)[0]
def timed(command, name, logs, env=None):
    logs.mkdir(parents=True,exist_ok=True)
    out=logs/(name+'.log'); tv=logs/(name+'.time-v.txt')
    start=time.monotonic()
    with out.open('wb') as f:
        r=subprocess.run(['/usr/bin/time','-v','-o',str(tv),'--',*command],stdout=f,stderr=subprocess.STDOUT,env=env)
    rss=None
    for line in tv.read_text().splitlines():
        if 'Maximum resident set size (kbytes)' in line: rss=int(line.split(':')[-1])*1024
    row=dict(command=command,exit_code=r.returncode,passed=r.returncode==0,wall_seconds=time.monotonic()-start,
             max_rss_bytes=rss,log=str(out),time_v=str(tv))
    return row
def baseline():
    require(sha(BASE/'BUILD_PLAN.json')==PLAN_SHA,'Baseline plan changed')
    require(sha(BASE/'BUILD_RESULT.json')==RESULT_SHA,'Baseline result changed')
    p,r=load(BASE/'BUILD_PLAN.json'),load(BASE/'BUILD_RESULT.json')
    require(r['passed'] is True and len(r['compile'])==158 and len(p['compile'])==158,
            'Baseline complete 158-TU link missing')
    require(r['plan_sha256']==PLAN_SHA and r['link']['passed'] and
            sha(r['candidate'])==r['candidate_sha256'],'Baseline actual plan/link/binary not intact')
    require(inventory(Path(p['source']))==p['source_manifest'],'Production baseline frozen source changed')
    verify_support(Path(p['source']),p['support_manifest'])
    require({x for x in p['full_link_command'] if x.endswith('.o')}==
            {row['object'] for row in r['compile']} and len({row['object'] for row in r['compile']})==158,
            'Baseline actual 158 objects differ from full production link')
    return p,r
def native_production_baseline():
    require(sha(NATIVE_BASE/'BUILD_PLAN.json')==NATIVE_PLAN_SHA,'Native recipe producer plan changed')
    require(sha(NATIVE_BASE/'BUILD_RESULT.json')==NATIVE_RESULT_SHA,'Native recipe producer result changed')
    return load(NATIVE_BASE/'BUILD_PLAN.json')
def native_recipes():
    path=NATIVE_BASE/'native_01/NATIVE_RESULT.json'
    require(sha(path)==NATIVE_RECIPE_SHA,'Baseline native recipes changed')
    result=load(path)
    require(len(result['tasks'])==48,'Expected complete 48 baseline native recipes')
    return result
def pin_native_external(external,p):
    for task in native_recipes()['tasks']:
        pin=task['dependency_file'];verify_pin(pin)
        source_root_prefix=p['source']+'/'
        for token in dep_tokens(pin['path']):
            if token.startswith(source_root_prefix):
                relative=token[len(source_root_prefix):]
                require(p['source_manifest'].get(relative)==sha(token),'Changed baseline native project dependency '+token)
            elif token.startswith(str(NATIVE_BASE)+'/native_01/') or '/build_attempt_03/native_01/' in token:
                # The generated wrapper/schema source is frozen separately;
                # no native object or validator from this directory is reused.
                continue
            else:
                require('/src/cpp/' not in token and '/src/Architecture/' not in token,
                        'Foreign native project dependency '+token)
                external.setdefault(token,sha(token))
def fresh_dependencies(plan,path):
    pins=[];prefix=plan['source']+'/'
    for token in dep_tokens(path):
        if token.startswith(prefix):
            relative=token[len(prefix):]
            expected=plan['source_manifest'].get(relative)
            if expected is None: expected=plan['support_manifest'].get(relative)
            if expected is None: expected=plan.get('generated_native_manifest',{}).get(relative)
            require(expected==sha(token),'Unfrozen fresh project/native dependency '+relative)
        else:
            require(plan['external_dependencies'].get(token)==sha(token),'Unpinned fresh dependency '+token)
        pins.append(file_pin(token))
    return dict(file=file_pin(path),members=pins)
def runtime_audit(executable,plan):
    completed=subprocess.run(['/usr/bin/ldd',str(executable)],capture_output=True,text=True,check=True)
    require('not found' not in completed.stdout,'Candidate has unresolved runtime library')
    actual=[];expected=plan['runtime_expected_by_resolved']
    for line in completed.stdout.splitlines():
        match=re.search(r'=>\s+(/\S+)',line) or re.match(r'\s*(/\S+)',line)
        if not match: continue
        pin=file_pin(match.group(1))
        require(expected.get(pin['resolved'])==pin['sha256'],'Unexpected/currently changed candidate runtime '+pin['path'])
        actual.append(pin)
    require(actual,'No actual candidate runtime libraries audited')
    return dict(command=['/usr/bin/ldd',str(executable)],output=completed.stdout,libraries=actual,verified_now=True)
def prior_current_native():
    for name,wanted in (('BUILD_PLAN.json',PRIOR_PLAN_SHA),('BUILD_RESULT.json',PRIOR_BUILD_SHA),
                        ('NATIVE_RESULT.json',PRIOR_NATIVE_SHA)):
        require(sha(PRIOR_NATIVE/name)==wanted,'Earlier current-candidate receipt changed '+name)
    p=load(PRIOR_NATIVE/'BUILD_PLAN.json');b=load(PRIOR_NATIVE/'BUILD_RESULT.json');n=load(PRIOR_NATIVE/'NATIVE_RESULT.json')
    require(b['passed'] and b['plan_sha256']==PRIOR_PLAN_SHA and b['link']['passed'] and
            b['candidate_sha256']==PRIOR_PRODUCER_SHA and sha(b['candidate'])==PRIOR_PRODUCER_SHA,
            'Earlier current production binary is not the verified producer')
    require(n['passed'] is False and n['plan_sha256']==PRIOR_PLAN_SHA and
            n['production_sha256']==PRIOR_PRODUCER_SHA and n['full_production_link'] is True and
            n['solver_run'] is False and len(n['tasks'])==50,'Earlier partial native scope changed')
    require(all(row.get('passed') for row in n['tasks'][:49]) and not n['tasks'][49]['passed'] and
            sum(bool(row.get('execute')) for row in n['tasks'][:49])==48 and
            n['tasks'][47]['execute_requested'] is False and 'execute' not in n['tasks'][47],
            'Expected 48 actual executions plus one build-only then compilation failure')
    for name,wanted in ANCESTOR_RECEIPTS.items():
        require(sha(ANCESTOR_NATIVE/name)==wanted,'Nested original current-candidate receipt changed '+name)
    ancestor_build=load(ANCESTOR_NATIVE/'BUILD_RESULT.json');ancestor=load(ANCESTOR_NATIVE/'NATIVE_RESULT.json')
    require(ancestor_build['passed'] and ancestor_build['candidate_sha256']==PRIOR_PRODUCER_SHA and
            sha(ancestor_build['candidate'])==PRIOR_PRODUCER_SHA and ancestor['passed'] is False and
            len(ancestor['tasks'])==47,'Nested original current-candidate scope changed')
    for index in range(46):
        adopted=n['tasks'][index];original=ancestor['tasks'][index]
        proof=adopted.get('acceptance_provenance',{})
        require(adopted.get('execution_reused_from_earlier_current_candidate') is True and
                proof.get('attempt')==str(ANCESTOR_NATIVE) and proof.get('record_index')==index+1 and
                proof.get('result_sha256')==ANCESTOR_RECEIPTS['NATIVE_RESULT.json'] and
                proof.get('source_aggregate_passed') is False,'Nested per-task adoption identity differs')
        for key in ('name','passed','execute_requested','test_sha256','args','compile','dependencies',
                    'object','link','executable','runtime','actual_args','execute'):
            require(adopted.get(key)==original.get(key),'Nested actual record was rewritten '+key)
        require(original['passed'] and all(original[phase]['passed'] and original[phase]['exit_code']==0
                for phase in ('compile','link','execute')),'Nested original actual execution did not pass')
    # The earlier adoption receipt carries literal file pins for its original
    # dependencies/artifacts; verify them now, not just their JSON labels.
    for field in ('earlier_current_candidate_pins','earlier_current_candidate_receipts',
                  'baseline_actual_project_pins','native_sidecars','tool_and_link_pins'):
        for pin in p.get(field,{}).values():verify_pin(pin)
    return p,b,n
def prior_native_acceptance(task,index,current_source,current_manifest,current_support,current_generated,
                            prior_plan,prior_build,prior_result,verified,pins):
    """Verify one historical EXECUTION, never accept the false aggregate."""
    row=prior_result['tasks'][index]
    require(row['name']==task['name'] and index<49 and row['passed'] and row['execute_requested']==task['execute'],
            'Wrong previous actual execution selected')
    phases=('compile','link','execute') if task['execute'] else ('compile','link')
    for phase in phases:
        require(row[phase]['passed'] and row[phase]['exit_code']==0,'Earlier task phase did not pass')
    if not task['execute']:require('execute' not in row,'Build-only record must not invent an execution')
    flags=row['compile']['command']
    require('-UNDEBUG' in flags and '-O0' in flags and
            max(i for i,x in enumerate(flags) if x=='-UNDEBUG')>
                max([-1,*[i for i,x in enumerate(flags) if x=='-DNDEBUG']]),
            'Earlier native assertions were not enabled')
    old_source=Path(row['compile']['command'][row['compile']['command'].index('-c')+1])
    require(sha(old_source)==row['test_sha256']==sha(task['source']),'Earlier/new test source bytes differ '+task['name'])
    require(row['args']==row.get('actual_args',row['args'])==task['args'] and
            (not task['execute'] or row['execute']['command']==[row['executable']['path'],*task['args']]),
            'Earlier/current test argv differ '+task['name'])
    def verify_once(pin):
        previous=verified.get(pin['path'])
        require(previous is None or previous['sha256']==pin['sha256'],'Conflicting earlier execution file identities')
        if previous is not None and 'resolved' in previous and 'resolved' in pin:
            require(previous['resolved']==pin['resolved'],'Conflicting earlier execution symlink identities')
        if previous is None:verify_pin(pin);verified[pin['path']]=pin
        pins[pin['path']]=pin
    verify_once(row['object']);verify_once(row['executable'])
    verify_once(row['dependencies']['file'])
    members=row['dependencies']['members']
    require(dep_tokens(row['dependencies']['file']['path'])==
            [os.path.normpath(pin['path']) for pin in members],'Earlier actual native .d differs from members')
    allowed_roots=(prior_plan['source'],str(ROOT/'build_03/source'))
    require(any(str(old_source).startswith(root+'/') for root in allowed_roots),
            'Earlier actual test source escaped the fixed 03/04 frozen roots')
    current_files={**current_manifest,**current_support,**current_generated}
    mapping=[]
    for pin in members:
        verify_once(pin)
        owner_root=next((root for root in allowed_roots if pin['path'].startswith(root+'/')),None)
        if owner_root is not None:
            relative=pin['path'][len(owner_root)+1:]
            require(current_files.get(relative)==pin['sha256'],
                    'Earlier/current dependency content differs '+relative)
            mapping.append(dict(previous=pin['path'],current=str(Path(current_source)/relative),sha256=pin['sha256']))
        else:
            require(prior_plan['external_dependencies'].get(pin['path'])==pin['sha256'],
                    'Earlier native dependency escaped frozen/external closure '+pin['path'])
    native_object=row['object']['path'];old_main=next(x['object'] for x in prior_build['compile'] if x['source']=='src/cpp/pea/main.cpp')
    renamed=str(Path(native_object).parent/'production_main_renamed.o')
    wanted={native_object,*[renamed if x['object']==old_main else x['object'] for x in prior_build['compile']]}
    actual=[x for x in row['link']['command'] if x.endswith('.o')]
    require(len(actual)==159 and set(actual)==wanted,'Earlier native link was not complete 158 producer objects')
    for pin in row['runtime']['libraries']:verify_once(pin)
    for phase in phases:
        for name in ('log','time_v'):
            value=file_pin(row[phase][name]);verify_once(value)
    for argument in task['args']:
        if Path(argument).is_file():verify_once(file_pin(argument))
    return dict(index=index+1,name=task['name'],dependency_mapping=mapping,
                previous_renamed_main=renamed,previous_test_source=file_pin(old_source),
                producer_sha256=PRIOR_PRODUCER_SHA,aggregate_was_failed=True,
                claim='INDIVIDUAL_ACTUAL_CURRENT_PRODUCER_EXECUTION' if task['execute'] else
                    'INDIVIDUAL_COMPILE_LINK_RUNTIME_BUILD_ONLY_NO_EXECUTION')
def baseline_dependencies(planned,actual,baseline_plan,external,verified,project_pins):
    """Resolve the OBJECT'S real .d/command ancestry, never guess build02 roots."""
    require(actual['source']==planned['source'] and actual['passed'],'Wrong baseline ordered TU')
    require(actual['object']==planned['object'],'Baseline planned/actual object mismatch')
    require(sha(actual['object'])==actual['object_sha256'],'Baseline object changed '+actual['source'])
    if not actual['reused']:
        dependency=actual['dependencies']['file']
        command=actual['command']
        require(command==planned['command'],'Fresh baseline command was not the planned command')
        recorded_members=actual['dependencies']['members']
    else:
        dependency=planned['prior_dep']
        command=planned['command']
        require(not planned['rebuild'] and actual['object']==planned['prior_object'] and
                actual['object_sha256']==planned['prior_object_sha256'],
                'Reused baseline object ancestry mismatch '+actual['source'])
        recorded_members=None
    verify_pin(dependency)
    tokens=dep_tokens(dependency['path'])
    if recorded_members is not None:
        require(tokens==[os.path.normpath(pin['path']) for pin in recorded_members],
                'Fresh baseline actual dependency list differs from audited members '+actual['source'])
        for pin in recorded_members:
            prior=verified.get(pin['path'])
            require(prior is None or prior['sha256']==pin['sha256'],'Conflicting baseline member hashes')
            if prior is None:verify_pin(pin);verified[pin['path']]=pin
    root=source_root(command);members=set()
    require(command[command.index('-c')+1]==root+'/'+actual['source'],
            'Actual dependency source root does not name the object TU '+actual['source'])
    for token in tokens:
        if token.startswith(root+'/'):
            relative=token[len(root)+1:]
            expected=baseline_plan['source_manifest'].get(relative)
            require(expected is not None,'Uncatalogued baseline actual project dependency '+token)
            pin=verified.get(token)
            if pin is None:
                pin=file_pin(token);verified[token]=pin
            require(pin['sha256']==expected,'Baseline actual ancestor source mismatch '+token)
            project_pins[token]=pin;members.add(relative)
        else:
            require('/src/cpp/' not in token and '/src/Architecture/' not in token,
                    'Foreign project dependency outside actual object root '+token)
            expected=baseline_plan['external_dependencies'].get(token)
            require(expected is not None,'Unpinned baseline external header '+token)
            pin=verified.get(token)
            if pin is None:pin=file_pin(token);verified[token]=pin
            require(pin['sha256']==expected,'Changed baseline external header '+token)
            external[token]=expected
    require(actual['source'] in members,'Object source absent from actual dependency closure')
    if actual['reused']:
        require(members==set(planned['project_dependencies']),
                'Reused baseline project closure does not match actual ancestor .d '+actual['source'])
    # A fresh build02 TU can legitimately introduce headers absent from its
    # older planning .d. Its ACTUAL recorded .d/members above are authoritative.
    return dependency,root,members,command
def prepare(attempt):
    require(os.uname().nodename.lower()!='comput','Local machine only')
    p,r=baseline();prior_plan,prior_build,prior_result=prior_current_native()
    dest=ROOT/attempt; ev=ART/attempt
    require(not dest.exists() and not ev.exists(),'New attempt required; preserve failed attempt')
    dest.mkdir(parents=True); ev.mkdir()
    src=dest/'source'; shutil.copytree(INPUT,src)
    manifest=inventory(src); require(manifest==inventory(INPUT),'Source changed while freezing')
    source_delta={key for key in set(manifest)|set(p['source_manifest'])
                  if manifest.get(key)!=p['source_manifest'].get(key)}
    expected_delta={
        'src/cpp/common/zhangCurrentEvidenceAuthoritySnapshot.hpp',
        'src/cpp/common/zhangNativeProductTargetAdmission.hpp',
        'src/cpp/common/zhangUs20E240Continuation.hpp',
        'src/cpp/pea/ppp_ambres.cpp',
        'src/cpp/pea/ppp_ambres_checkpoint.inc',
        'src/cpp/pea/zhangPppAr.cpp',
        'src/tests/unit/current_evidence_authority_snapshot.cpp',
        'src/tests/unit/gec_ar_runtime_checkpoint.cpp',
        'src/tests/unit/native_product_target_admission.cpp',
        'src/tests/unit/us20_e240_continuation_policy.cpp'}
    require(source_delta==expected_delta,'Authority repair source scope differs from the reviewed ten-file delta')
    save(ev/'SOURCE_DELTA.json',dict(schema='AUTHORITY_REPAIR_SOURCE_DELTA_V1',passed=True,
        baseline_plan=file_pin(BASE/'BUILD_PLAN.json'),changed_files=[dict(path=key,
            baseline_sha256=p['source_manifest'].get(key),new_sha256=manifest.get(key)) for key in sorted(source_delta)],
        numerical_solver_sources_unchanged=False,full_covariance_and_safety_checks_retained=True))
    support_manifest=support_inventory(src)
    require({Path(k).name for k in support_manifest}==SUPPORT_FILES,'Exactly three codec support inputs required')
    require(support_manifest==support_inventory(INPUT),'Codec support changed while freezing')
    support_provenance=verify_support(src,support_manifest)
    prior_changes=[k for k in sorted(set(manifest)|set(prior_plan['source_manifest']))
                   if manifest.get(k)!=prior_plan['source_manifest'].get(k)]
    reuse_current_executions=False  # All native acceptance must use this new producer.
    build=dest/'candidate'; remap={}; recipes=[]; external={};verified={};project_pins={}
    for planned,actual in zip(p['compile'],r['compile']):
        pin,root,members,old=baseline_dependencies(planned,actual,p,external,verified,project_pins)
        changed=sorted(k for k in members if manifest.get(k)!=p['source_manifest'].get(k))
        old_obj=actual['object']
        obj=str(build/'objects'/Path(actual['source']).with_suffix('.o')) if changed else actual['object']
        command_root=source_root(old)
        cmd=[x.replace(command_root,str(src)) for x in old]
        for flag in ('-o','-MF'):
            cmd[cmd.index(flag)+1]=obj+('.d' if flag=='-MF' else '')
        # -MT is diagnostic only; keep it local as well.
        if '-MT' in cmd: cmd[cmd.index('-MT')+1]=obj
        recipes.append(dict(source=actual['source'],object=obj,command=cmd if changed else old,
                            rebuild=bool(changed),changed_dependencies=changed,project_dependencies=sorted(members),
                            baseline_actual_source_root=root,
                            prior_object=actual['object'],prior_object_sha256=actual['object_sha256'],
                            prior_dep=pin))
        remap[old_obj]=obj
    # Native-only system headers must be frozen before any native compilation,
    # too. Old native dependencies are recipes/provenance, never accepted test
    # results or reusable native objects.
    native_p=native_production_baseline()
    pin_native_external(external,native_p)
    generated_manifest={};native_tasks=[];native_sidecars={}
    generated_dir=src/'native_generated';generated_dir.mkdir()
    old_native=native_recipes()
    for task in old_native['tasks']:
        name=task['name'];arguments=task.get('execute',{}).get('command',[])[1:]
        for argument in arguments:
            if Path(argument).is_file(): native_sidecars[argument]=file_pin(argument)
        if name=='ratio_only_skip_nis_native':
            original=Path(task['compile']['command'][task['compile']['command'].index('-c')+1])
            original_pin=file_pin(original);text=original.read_text()
            require(str(native_p['source']) in text,'Generated ratio wrapper source prefix unrecognised')
            if reuse_current_executions:
                previous=prior_result['tasks'][45]
                require(previous['name']==name,'Earlier ratio execution index changed')
                previous_source=Path(previous['compile']['command'][previous['compile']['command'].index('-c')+1])
                require(str(previous_source).startswith(str(ROOT/'build_03/source')+'/'),
                        'Earlier ratio wrapper must be the exact original build03 frozen source')
                require(sha(previous_source)==previous['test_sha256'],'Earlier ratio wrapper changed')
                # Exact bytes of the actually executed wrapper, including its
                # old frozen include. EVERY included old file is separately
                # compared with the current frozen relative file below.
                text=previous_source.read_text()
            else:
                text=text.replace(str(native_p['source']),str(src))
            target=generated_dir/original.name;target.write_text(text)
            generated_manifest[target.relative_to(src).as_posix()]=sha(target)
            native_tasks.append(dict(name=name,source=str(target),execute=True,args=arguments,
                                     original_generated_source=original_pin))
        else:
            relative='src/tests/unit/'+name+'.cpp';target=src/relative
            require(relative in manifest and target.is_file(),'Missing frozen baseline native source '+name)
            changed=manifest[relative]!=native_p['source_manifest'].get(relative)
            native_tasks.append(dict(name=name,source=str(target),execute=False if name=='us20_e240_nine_section_audit' else task['execute_requested'] or changed,args=arguments))
    included={task['name'] for task in native_tasks}
    for relative,digest in manifest.items():
        if (relative.startswith('src/tests/unit/') and relative.endswith('.cpp') and
            native_p['source_manifest'].get(relative)!=digest and Path(relative).stem not in included):
            name=Path(relative).stem
            native_tasks.append(dict(name=name,source=str(src/relative),execute=True,args=[]));included.add(name)
    tool_pins={}
    for path,wanted in p['tool_and_link_pins'].items():
        verify_pin(wanted);tool_pins[path]=file_pin(path)
    for path in (str(Path(__file__).resolve()),sys.executable,'/usr/bin/objcopy','/usr/bin/ldd','/usr/bin/time'):
        tool_pins[path]=file_pin(path)
    runtime_expected={}
    for pin in r['runtime']['libraries']:
        verify_pin(pin);current=file_pin(pin['path']);tool_pins[pin['path']]=current
        runtime_expected[current['resolved']]=current['sha256']
    for path,pin in native_sidecars.items():verify_pin(pin)
    link=[remap.get(x,x) for x in p['full_link_command']]
    link[link.index('-o')+1]=str(build/'pea')
    require(sum(x.endswith('.o') for x in link)==158,'Full 158 TU link required')
    require({x for x in link if x.endswith('.o')}=={row['object'] for row in recipes} and
            len({row['object'] for row in recipes})==158,'Full link objects do not equal exact 158-TU closure')
    source_changes=[dict(path=k,baseline_sha256=p['source_manifest'].get(k),new_sha256=manifest.get(k))
                    for k in sorted(set(manifest)|set(p['source_manifest']))
                    if manifest.get(k)!=p['source_manifest'].get(k)]
    test_only_delta=all(change['path'].startswith('src/tests/unit/') and change['path'].endswith('.cpp')
                        for change in source_changes)
    if test_only_delta:require(not any(row['rebuild'] for row in recipes),'Test-only delta forced production rebuild')
    require(len(native_tasks)==53,'Expected all 51 baseline cases plus both authority repair regressions')
    accepted={};prior_pins={}
    if reuse_current_executions:
        for index,task in enumerate(native_tasks[:49]):
            accepted[task['name']]=prior_native_acceptance(task,index,str(src),manifest,support_manifest,
                generated_manifest,prior_plan,prior_build,prior_result,verified,prior_pins)
    plan=dict(schema='DATUM_UNKNOWN_ISOLATED_BUILD_V4',source=str(src),build=str(build),evidence=str(ev),
              source_manifest=manifest,support_manifest=support_manifest,support_provenance=support_provenance,
              generated_native_manifest=generated_manifest,native_tasks=native_tasks,native_sidecars=native_sidecars,
              native_recipe_pin=file_pin(NATIVE_BASE/'native_01/NATIVE_RESULT.json'),
              native_producer_plan_pin=file_pin(NATIVE_BASE/'BUILD_PLAN.json'),
              baseline_actual_project_pins=project_pins,source_changes=source_changes,test_only_source_delta=test_only_delta,
              earlier_current_candidate_acceptance=accepted,earlier_current_candidate_pins=prior_pins,
              earlier_current_candidate_receipts={name:file_pin(PRIOR_NATIVE/name) for name in
                  ('BUILD_PLAN.json','BUILD_RESULT.json','NATIVE_RESULT.json')},
              nested_current_candidate_receipts={name:file_pin(ANCESTOR_NATIVE/name) for name in ANCESTOR_RECEIPTS},
              source_changes_since_prior_current_native=prior_changes,
              external_dependencies=external,tool_and_link_pins=tool_pins,runtime_expected_by_resolved=runtime_expected,
              compile=recipes,full_link_command=link,
              runtime_lib_directory=p['runtime_lib_directory'],baseline_plan_sha256=PLAN_SHA,baseline_result_sha256=RESULT_SHA,
              serial=True,solver_run=False,server_access=False,ratio=2,joint_nis=False,epoch_independent_ar=True)
    save(ev/'BUILD_PLAN.json',plan)
    print(json.dumps(dict(prepared=True,plan_sha256=sha(ev/'BUILD_PLAN.json'),fresh=sum(x['rebuild'] for x in recipes),
        reused=sum(not x['rebuild'] for x in recipes),individually_verified_prior_current_records=len(accepted),
        individually_verified_prior_current_executions=sum(t['execute'] for t in native_tasks if t['name'] in accepted),
        individually_verified_prior_current_build_only=sum(not t['execute'] for t in native_tasks if t['name'] in accepted),
        source_changes=[x['path'] for x in source_changes])),flush=True)
def verify_frozen(plan):
    require(inventory(Path(plan['source']))==plan['source_manifest'],'Frozen source changed')
    verify_support(Path(plan['source']),plan['support_manifest'])
    for relative,wanted in plan.get('generated_native_manifest',{}).items():
        require(sha(Path(plan['source'])/relative)==wanted,'Frozen generated native wrapper changed '+relative)
    for path,wanted in plan['external_dependencies'].items(): require(sha(path)==wanted,'External dependency changed '+path)
    for pin in plan['tool_and_link_pins'].values():verify_pin(pin)
    for pin in plan['native_sidecars'].values():verify_pin(pin)
    for pin in plan.get('baseline_actual_project_pins',{}).values():verify_pin(pin)
    for row in plan['compile']:verify_pin(row['prior_dep'])
    verify_pin(plan['native_recipe_pin'])
    verify_pin(plan['native_producer_plan_pin'])
    for pin in plan.get('earlier_current_candidate_pins',{}).values():verify_pin(pin)
    for pin in plan.get('earlier_current_candidate_receipts',{}).values():verify_pin(pin)
    for pin in plan.get('nested_current_candidate_receipts',{}).values():verify_pin(pin)
    for task in plan['native_tasks']:
        if 'original_generated_source' in task:verify_pin(task['original_generated_source'])
def build(attempt):
    ev=ART/attempt; plan=load(ev/'BUILD_PLAN.json'); verify_frozen(plan)
    dest=Path(plan['build']); require(not dest.exists(),'Build is single-use')
    dest.mkdir(); (dest/'lib').symlink_to(plan['runtime_lib_directory'],target_is_directory=True)
    result=dict(passed=False,compile=[],plan_sha256=sha(ev/'BUILD_PLAN.json'),solver_run=False)
    try:
        for i,row in enumerate(plan['compile']):
            if row['rebuild']:
                Path(row['object']).parent.mkdir(parents=True,exist_ok=True)
                print('COMPILE',row['source'],flush=True)
                record=timed(row['command'],f'compile_{i:03d}_{Path(row["source"]).stem}',ev/'logs')
                record.update(source=row['source'],object=row['object'],reused=False)
                result['compile'].append(record); save(ev/'progress.json',result)
                require(record['passed'],'Compile failed '+row['source'])
                result['compile'][-1]['dependencies']=fresh_dependencies(plan,row['object']+'.d')
            else:
                require(sha(row['object'])==row['prior_object_sha256'],'Reused object changed')
                result['compile'].append(dict(source=row['source'],object=row['object'],reused=True,passed=True))
            result['compile'][-1]['object_sha256']=sha(row['object'])
        result['link']=timed(plan['full_link_command'],'full_158_object_link',ev/'logs')
        require(result['link']['passed'],'Complete link failed')
        result['runtime']=runtime_audit(dest/'pea',plan)
        verify_frozen(plan)
        result.update(passed=True,candidate=str(dest/'pea'),candidate_sha256=sha(dest/'pea'))
    except BaseException as e:
        result['error']=repr(e); raise
    finally: save(ev/'BUILD_RESULT.json',result)
    print(json.dumps(dict(passed=True,candidate=result['candidate'],sha256=result['candidate_sha256'])),flush=True)
def compile_probe(name):
    p,r=baseline(); old=p['compile'][0]['command']; root=source_root(old)
    directory=ROOT/'header_probes'; directory.mkdir(parents=True,exist_ok=True)
    obj=str(directory/(name+'.o'))
    cmd=[x.replace(root,str(INPUT)) for x in old]
    for flag,value in [('-c',str(INPUT/'src/tests/unit'/(name+'.cpp'))),('-o',obj),('-MF',obj+'.d'),('-MT',obj)]: cmd[cmd.index(flag)+1]=value
    cmd.insert(cmd.index('-c'),'-UNDEBUG')
    answer=timed(cmd,name,ART/'header_probes'); save(ART/'header_probes'/(name+'.json'),answer)
    print(json.dumps(answer,indent=2)); require(answer['passed'],'Header probe failed')
def native(attempt):
    ev=ART/attempt; plan=load(ev/'BUILD_PLAN.json'); built=load(ev/'BUILD_RESULT.json')
    require(built['passed'] and sha(built['candidate'])==built['candidate_sha256'],'Production link required')
    verify_frozen(plan)
    for record in built['compile']:
        require(sha(record['object'])==record['object_sha256'],'Production object changed before native link '+record['source'])
    prior_plan,prior_build,prior_result=prior_current_native()
    accepted=plan.get('earlier_current_candidate_acceptance',{})
    if accepted:
        require(built['candidate_sha256']==prior_result['production_sha256']==PRIOR_PRODUCER_SHA,
                'Cannot reuse executions from a different producer')
        require({row['source']:row['object_sha256'] for row in built['compile']}==
                {row['source']:row['object_sha256'] for row in prior_build['compile']},
                'Current and earlier actual production object SHA sets differ')
    output=Path(plan['build']).parent/'native_01'; require(not output.exists(),'New native attempt required')
    output.mkdir(); (output/'lib').symlink_to(plan['runtime_lib_directory'],target_is_directory=True)
    # The old baseline03 aggregate supplies recipes only. Separately, the
    # current build04 producer's first 49 actual records have an explicit
    # per-record acceptance path (48 executions plus one build-only); the
    # original build03 nested acceptance chain is also checked literally.
    tasks=plan['native_tasks']
    main=next(r['object'] for r in plan['compile'] if r['source']=='src/cpp/pea/main.cpp')
    renamed=output/'production_main_renamed.o'
    rename=timed(['/usr/bin/objcopy','--redefine-sym','main=zhang_native_production_main',main,str(renamed)],'rename_main',ev/'native_logs')
    require(rename['passed'],'Native main renaming failed')
    result=dict(schema='DATUM_UNKNOWN_NATIVE_AGGREGATE_V4',passed=False,tasks=[],plan_sha256=sha(ev/'BUILD_PLAN.json'),production_sha256=built['candidate_sha256'],
                full_production_link=True,whole_network_real_replay=False,solver_run=False)
    result['cross_attempt_acceptance']=dict(source_attempt=str(PRIOR_NATIVE),source_result_sha256=PRIOR_NATIVE_SHA,
        source_aggregate_passed=False,only_individual_actual_pass_records_reusable=True,
        producer_sha256=PRIOR_PRODUCER_SHA,
        expected_reused_executions=sum(task['execute'] for task in tasks if task['name'] in accepted),
        expected_reused_build_only=sum(not task['execute'] for task in tasks if task['name'] in accepted),
        nested_ancestor_attempt=str(ANCESTOR_NATIVE),nested_ancestor_receipt_sha256=ANCESTOR_RECEIPTS['NATIVE_RESULT.json'])
    source=Path(plan['source']); recipe=plan['compile'][0]['command']; recipe_root=source_root(recipe)
    # The private codec header contains declarations only, not validator code.
    support=source/'test_support'
    env=os.environ.copy(); env.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
    try:
        for task in tasks:
            name=task['name'];test=Path(task['source']);execute=task['execute']
            print('NATIVE',name,flush=True)
            if name in accepted:
                evidence=accepted[name];previous=prior_result['tasks'][evidence['index']-1]
                require(previous['name']==name and previous['passed'] and
                        sha(test)==previous['test_sha256'] and
                        task['args']==previous.get('actual_args',previous['args']),
                        'Earlier actual execution no longer matches current task '+name)
                require(sha(evidence['previous_renamed_main'])==sha(renamed),
                        'Earlier renamed main object differs from current actual producer')
                actual_runtime=runtime_audit(previous['executable']['path'],plan)
                require(sorted((x['resolved'],x['sha256']) for x in actual_runtime['libraries'])==
                        sorted((x['resolved'],x['sha256']) for x in previous['runtime']['libraries']),
                        'Earlier execution runtime closure changed '+name)
                row=dict(previous)
                row['acceptance_reused_from_earlier_current_candidate']=True
                row['execution_reused_from_earlier_current_candidate']=bool(execute)
                row['build_only_reused_from_earlier_current_candidate']=not execute
                row['current_source']=file_pin(test)
                row['acceptance_provenance']=dict(attempt=str(PRIOR_NATIVE),result_sha256=PRIOR_NATIVE_SHA,
                    record_index=evidence['index'],source_aggregate_passed=False,
                    exact_source_sha256_match=True,actual_dependency_closure_verified=True,
                    identical_producer_and_158_object_sha_set=True,
                    object_executable_argv_and_sidecars_verified=True,
                    current_runtime=actual_runtime)
                result['tasks'].append(row);save(ev/'native_progress.json',result)
                continue
            obj=output/(name+'.o'); command=[x.replace(recipe_root,str(source)) for x in recipe]
            for flag,value in (('-c',str(test)),('-o',str(obj)),('-MF',str(obj)+'.d'),('-MT',str(obj))): command[command.index(flag)+1]=value
            command[command.index('-c'):command.index('-c')]=['-UNDEBUG','-O0','-I'+str(support)]
            if name=='historical_evidence_archive':command.insert(command.index('-c'),'-DZHANG_ARCHIVE_PRODUCTION_LINK=1')
            if name=='us20_e240_continuation_policy':command.insert(command.index('-c'),'-DZHANG_E240_POLICY_PRODUCTION_LINK=1')
            row=dict(name=name,passed=False,execute_requested=execute,test_sha256=sha(test),args=task['args'])
            result['tasks'].append(row)
            row['compile']=timed(command,name+'_compile',ev/'native_logs',env)
            require(row['compile']['passed'],'Native compile failed '+name)
            row['dependencies']=fresh_dependencies(plan,str(obj)+'.d')
            row['object']=file_pin(obj)
            link=[plan['full_link_command'][0],str(obj),*plan['full_link_command'][1:]]
            link=[str(renamed) if x==main else x for x in link]
            link[link.index('-o')+1]=str(output/name)
            row['link']=timed(link,name+'_link',ev/'native_logs',env)
            require(row['link']['passed'],'Native link failed '+name)
            row['executable']=file_pin(output/name)
            row['runtime']=runtime_audit(output/name,plan)
            if execute:
                arguments=task['args']
                if name=='historical_evidence_archive':
                    evidence=output/(name+'_evidence')
                    require(not evidence.exists(),'Archive test requires new independent evidence directory')
                    arguments=[str(evidence)]
                row['actual_args']=arguments
                row['execute']=timed([str(output/name),*arguments],name+'_run',ev/'native_logs',env)
                require(row['execute']['passed'],'Native regression failed '+name)
            row['passed']=True; save(ev/'native_progress.json',result)
        verify_frozen(plan)
        reuse_count=sum(bool(x.get('execution_reused_from_earlier_current_candidate')) for x in result['tasks'])
        new_count=sum(bool(x.get('execute')) for x in result['tasks'] if not x.get('execution_reused_from_earlier_current_candidate'))
        require(len(result['tasks'])==53 and all(x['passed'] for x in result['tasks']),
                'Complete independent 53-case aggregate required')
        result.update(passed=True,run_cases=sum(x['execute'] for x in tasks),
            actual_executions_in_this_attempt=new_count,accepted_earlier_current_candidate_executions=reuse_count,
            build_only_cases=sum(not x['execute'] for x in tasks),selected_cases=53,
            accepted_earlier_build_only_records=sum(bool(x.get('build_only_reused_from_earlier_current_candidate')) for x in result['tasks']),
            acceptance_description='48 earlier individually verified current-candidate executions plus 2 new executions plus 1 earlier build-only'
                if reuse_count==48 and new_count==2 else 'Explicitly counted earlier and new individual executions; no failed aggregate promoted')
    except BaseException as e:
        result['error']=repr(e); raise
    finally: save(ev/'NATIVE_RESULT.json',result)
    print(json.dumps(dict(passed=True,run_cases=result['run_cases'],build_only_cases=result['build_only_cases'],
        prior_current_candidate_executions=result['accepted_earlier_current_candidate_executions'],
        actual_executions_in_this_attempt=result['actual_executions_in_this_attempt'])),flush=True)
if __name__=='__main__':
    a=argparse.ArgumentParser(); a.add_argument('mode',choices=['prepare','build','probe','native']); a.add_argument('name'); args=a.parse_args()
    {'prepare':prepare,'build':build,'probe':compile_probe,'native':native}[args.mode](args.name)
