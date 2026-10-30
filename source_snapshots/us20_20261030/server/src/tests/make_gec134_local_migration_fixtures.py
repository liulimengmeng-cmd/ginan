"""Independent pure-policy fixtures from SHA-locked real source canonical bytes.

No production policy/builder is imported, no estimator is run, and no source
archive is rewritten. Output is a new fixture file, never an old evidence file.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path

BASE = '/home/rx/GINAN/multi_gnss_float_profile_20261022'
SOURCE = BASE + '/server134_fixedcoord_float1h_offline_v3_20261023_01/checkpoints/E29-20240717010000-e121-r32f626de9463'
LOCAL = BASE + '/server134_fixedcoord_local_memory_20261024_01'
FLAG = 'ZHANG_GEC134_LOCAL_ONE_EPOCH_12_6'
POLICY = 'GEC134_FLOAT_E121_TO_AR_LOCAL_12_6_ONE_EPOCH_V1'
SOURCE_HASH = '6ec5bbb76c9c0cc977b7c2bcc5197a44b3ab70ff239cb6a418ac5fd9a085b21a'
SOURCE_BUNDLE = '1310d0f6543e552ac318208b79491139d93886de87d8331ee14cb14dc3bf283b'


def scalar(key):
    return key in {'command_argument_count', 'config_file_count', 'effective_capture_epoch_count'} or (
        key.startswith('config_file_') and key.endswith('_size'))


def decode(text):
    data = text.encode('utf-8')
    prefix = b'E29_CANONICAL_CONFIG_V1\n'
    assert data.startswith(prefix)
    cursor, result = len(prefix), {}
    while cursor < len(data):
        end = data.index(b'\n', cursor)
        key, length = data[cursor:end].split(b'=', 1)
        key, size = key.decode('ascii'), int(length)
        cursor = end+1
        if scalar(key):
            value = length.decode('ascii')
        else:
            value = data[cursor:cursor+size].decode('utf-8')
            assert data[cursor+size:cursor+size+1] == b'\n'
            cursor += size+1
        assert key not in result
        result[key] = value
    return result


def encode(records):
    parts = ['E29_CANONICAL_CONFIG_V1\n']
    for key, value in records.items():
        value = str(value)
        parts.append(f'{key}={value}\n' if scalar(key) else f'{key}={len(value.encode())}\n{value}\n')
    return ''.join(parts)


def overlay(run, end='2024-07-17 01:00:30', maximum=122):
    return {'processing_options': {'epoch_control': {'end_epoch': end, 'max_epochs': maximum},
        'gnss_general': {'zhang_pppar': {'checkpoint_capture_epochs': [end],
            'checkpoint_output_directory': run+'/checkpoints', 'checkpoint_restore_path': SOURCE,
            'product_filename': run+'/outputs/zhang_internal_products.csv',
            'product_covariance_filename': run+'/outputs/zhang_internal_product_covariance.csv'}}},
        'outputs': {'outputs_root': run+'/outputs', 'metadata': {
            'config_description': 'INDEPENDENT_POLICY_FIXTURE_NOT_A_SOLVER_RUN'}}}


def make_target(source, node, local=False):
    target = copy.deepcopy(source)
    del target['environment_ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT']
    for key, value in {'ZHANG_MULTI_GNSS_AR_CHECKPOINT': '1', 'ZHANG_GEC134_FLOAT_TO_AR_MIGRATION': '1',
        'ZHANG_FLOAT_MATRIX_THREADS': '6' if local else '14', 'ZHANG_FLOAT_MEMORY_MODE': 'parallel',
        'ZHANG_FLOAT_PURE_OBS_RANK_DIAGNOSTIC': 'off', 'ZHANG_FLOAT_PROFILE': '1'}.items():
        target['environment_'+key] = value
    target['environment_ZHANG_R51_AR_START_GPST_SECONDS'] = '1405213230'
    if local:
        target.update({'environment_'+FLAG: '1', 'gec134_local_migration_policy': POLICY,
            'environment_OMP_NUM_THREADS': '12', 'environment_ZHANG_R51_FLOAT_THREADS': '12',
            'float_available_cpus': '12', 'float_physical_cores': '6'})
    raw = json.dumps(node, ensure_ascii=False, indent=2)+'\n'
    run = node['outputs']['outputs_root'].removesuffix('/outputs')
    path = run+'/config/99_gec134_float_to_ar_overlay.yaml'
    target.update({'config_file_count': '11', 'config_file_10_path': path,
        'config_file_10_size': str(len(raw.encode())), 'config_file_10_bytes': raw,
        'config_file_10_sha256': hashlib.sha256(raw.encode()).hexdigest()})
    args = [source[f'command_argument_{i}'] for i in range(1,150)]
    args.insert(10, path)
    for i in range(1,150): del target[f'command_argument_{i}']
    target['command_argument_count'] = '150'
    target.update({f'command_argument_{i+1}': arg for i,arg in enumerate(args)})
    epoch = node['processing_options']['epoch_control']
    ar = node['processing_options']['gnss_general']['zhang_pppar']
    target.update({'effective_end_epoch': epoch['end_epoch'].replace(' ','T'),
        'effective_max_epochs': str(epoch['max_epochs']), 'effective_outputs_root': node['outputs']['outputs_root'],
        'effective_checkpoint_output_directory': ar['checkpoint_output_directory'],
        'effective_checkpoint_restore_path': ar['checkpoint_restore_path'],
        'effective_product_filename': ar['product_filename'],
        'effective_product_covariance_filename': ar['product_covariance_filename']})
    for i in range(int(target['effective_capture_epoch_count'])): del target[f'effective_capture_epoch_{i}']
    target['effective_capture_epoch_count'] = str(len(ar['checkpoint_capture_epochs']))
    target.update({f'effective_capture_epoch_{i}': value for i,value in enumerate(ar['checkpoint_capture_epochs'])})
    return target


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source_manifest', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    manifest = json.loads(args.source_manifest.read_text(encoding='utf-8'))
    text = manifest['config_text']
    assert manifest['config_sha256'] == SOURCE_HASH == hashlib.sha256(text.encode()).hexdigest()
    source = decode(text)
    assert encode(source) == text
    # Values independently checked against the actual locked archive, NOT
    # guessed from the number of logical processors or copied from new policy.
    expected_cpu = {'float_available_cpus':'28', 'float_physical_cores':'28',
        'environment_OMP_NUM_THREADS':'28', 'environment_ZHANG_R51_FLOAT_THREADS':'28'}
    assert all(source[key] == value for key,value in expected_cpu.items())
    assert source['working_directory'] == '/mnt/d/GINAN_R20/inputData'
    assert source['config_file_count'] == '10' and source['command_argument_count'] == '149'
    env = {'ZHANG_GEC134_FLOAT_TO_AR_MIGRATION':'1', 'ZHANG_MULTI_GNSS_AR_CHECKPOINT':'1',
        'ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT':'UNSET', 'ZHANG_R51_REUSE_FLOAT_20240717':'UNSET',
        'ZHANG_R51_REUSE_RATIO_20240717':'UNSET', 'ZHANG_R51_ENABLE':'1', 'ZHANG_R51_RATIO_ONLY':'1',
        'ZHANG_R51_RATIO_THRESHOLD':'2', 'ZHANG_R51_EPOCH_ONLY_AR':'1',
        'ZHANG_R51_AR_START_GPST_SECONDS':'1405213230', 'ZHANG_R51_FLOAT_THREADS':'28',
        'OMP_NUM_THREADS':'28', 'OPENBLAS_NUM_THREADS':'1', 'MKL_NUM_THREADS':'1',
        'ZHANG_FLOAT_MATRIX_THREADS':'14', 'ZHANG_FLOAT_MEMORY_MODE':'parallel',
        'ZHANG_FLOAT_PURE_OBS_RANK_DIAGNOSTIC':'off', 'ZHANG_FLOAT_PROFILE':'1',
        'ZHANG_P0_RESOURCE_PROBE':'1', 'ZHANG_R49_FUSION_WEIGHT':'0.10',
        'ZHANG_R51_SINGLE_BASIS_HANDOFF':'1', 'ZHANG_R51_PRODUCT_REPAIR':'1',
        'ZHANG_R51_INTERN_HANDOFF_PROOFS':'1', 'ZHANG_E29_RESTORE_BUNDLE': SOURCE}
    local_env = {FLAG:'1', 'ZHANG_R51_FLOAT_THREADS':'12', 'OMP_NUM_THREADS':'12', 'ZHANG_FLOAT_MATRIX_THREADS':'6'}
    normal = make_target(source, overlay(BASE+'/server134_fixedcoord_policy_fixture'))
    local = make_target(source, overlay(LOCAL), True)
    cases = []
    def add(name, expected=False, records=None, profile='server', env_overrides=None, **kwargs):
        changes = dict(local_env) if profile == 'local' else {}
        changes.update(env_overrides or {})
        cases.append(dict(name=name, expected=expected, profile=profile,
            target=encode((local if profile == 'local' else normal) if records is None else records),
            env_overrides=changes, **kwargs))
    add('server_unchanged_first_epoch', True)
    add('server_unchanged_restore_only', True, env_overrides={'ZHANG_GEC134_RESTORE_ONLY':'1'})
    add('server_unchanged_full_hour', True,
        make_target(source, overlay(BASE+'/server134_fixedcoord_policy_fixture','2024-07-17 02:00:00',241)))
    add('local_exact_one_epoch', True, profile='local')
    add('local_exact_restore_only', True, profile='local', env_overrides={'ZHANG_GEC134_RESTORE_ONLY':'1'})
    for profile in ('server', 'local'):
        reference = local if profile == 'local' else normal
        for key in ('source_bundle_sha','source_binary_sha','source_config_sha','source_input_sha','target_input_sha'):
            add(profile+'_wrong_'+key, profile=profile, **{key:'0'*64})
        for key,value in [('ZHANG_R51_RATIO_ONLY','0'),('ZHANG_R51_RATIO_THRESHOLD','3'),
            ('ZHANG_R51_EPOCH_ONLY_AR','0'),('ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT','1'),
            ('ZHANG_R51_AR_START_GPST_SECONDS','1405213200'),('ZHANG_FLOAT_PURE_OBS_RANK_DIAGNOSTIC','svd'),
            ('ZHANG_R51_JOINT_NIS','1'),('ZHANG_R51_NIS_ENABLED','0'),('ZHANG_NEW_UNKNOWN_OVERRIDE','1')]:
            add(profile+'_bad_environment_'+key, profile=profile, env_overrides={key:value})
        for key,value in [('config_file_9_bytes','different coordinates'),('config_file_9_sha256','0'*64),
            ('config_file_0_path','/different/boxwing.yaml'),('command_argument_15','/different/obs.rnx'),
            ('effective_start_epoch','2024-07-17T01:00:30'),('effective_checkpoint_runtime_id','OTHER_OWNER'),
            ('working_directory','/different/inputData'),('unknown_canonical_record','1')]:
            changed=reference.copy();changed[key]=value
            add(profile+'_bad_canonical_'+key, records=changed, profile=profile)
        for key,value in [('epoch_interval',60),('wait_next_epoch',30.05),('start_epoch','2024-07-17 01:00:30')]:
            node=overlay(LOCAL if profile=='local' else BASE+'/server134_fixedcoord_policy_fixture')
            node['processing_options']['epoch_control'][key]=value
            add(profile+'_forbidden_overlay_'+key, records=make_target(source,node,profile=='local'), profile=profile)
        node=overlay(LOCAL if profile=='local' else BASE+'/server134_fixedcoord_policy_fixture')
        node['receiver_options']={'MAL2':{'apriori_position':[1,2,3]}}
        add(profile+'_changed_station_model', records=make_target(source,node,profile=='local'), profile=profile)
        add(profile+'_changed_source_bytes',profile=profile,source_config_text=text+'\n')
        add(profile+'_missing_ratio_env',profile=profile,remove_env=['ZHANG_R51_RATIO_ONLY'])
    for key,value in [('OMP_NUM_THREADS','28'),('ZHANG_R51_FLOAT_THREADS','all'),('ZHANG_FLOAT_MATRIX_THREADS','14')]:
        add('local_wrong_thread_'+key,profile='local',env_overrides={key:value})
    for key,value in [('environment_OMP_NUM_THREADS','28'),('environment_ZHANG_R51_FLOAT_THREADS','28'),
        ('float_available_cpus','28'),('float_available_cpus','6'),('float_physical_cores','12'),
        ('float_physical_cores','28'),('gec134_local_migration_policy','GEC134_FLOAT_E121_TO_AR_V1'),
        ('environment_'+FLAG,'UNSET')]:
        changed=local.copy();changed[key]=value
        add('local_wrong_actual_'+key+'_'+value,records=changed,profile='local')
    for key in ('float_available_cpus','float_physical_cores','gec134_local_migration_policy','environment_'+FLAG):
        changed=local.copy();del changed[key]
        add('local_missing_actual_'+key,records=changed,profile='local')
    for value in ('0','true','UNSET'):
        add('local_bad_explicit_flag_'+value,profile='local',env_overrides={FLAG:value})
    add('local_flag_missing_from_environment',profile='local',remove_env=[FLAG])
    for end, maximum in [('2024-07-17 01:01:00',123),('2024-07-17 02:00:00',241),
        ('2024-07-17 01:00:00',121),('2024-07-17 01:00:30',123)]:
        add('local_not_exact_one_epoch_'+end+'_'+str(maximum),profile='local',
            records=make_target(source,overlay(LOCAL,end,maximum),True))
    for captures in ([], ['2024-07-17 01:00:30']*2,
        ['2024-07-17 01:00:00','2024-07-17 01:00:30'],['2024-07-17 01:01:00']):
        node=overlay(LOCAL);node['processing_options']['gnss_general']['zhang_pppar']['checkpoint_capture_epochs']=captures
        add('local_bad_captures_'+repr(captures),profile='local',records=make_target(source,node,True))
    for root in (LOCAL+'_02',BASE+'/server134_fixedcoord_other_local',BASE+'/../outside',
        BASE+'/server134_fixedcoord_float1h_offline_v3_20261023_01'):
        add('local_forbidden_root_'+root,profile='local',records=make_target(source,overlay(root),True))
    changed=normal.copy();changed['float_available_cpus']='12'
    add('server_local_cpu_without_flag_rejected',records=changed)
    changed=normal.copy();changed['environment_'+FLAG]='1';changed['gec134_local_migration_policy']=POLICY
    add('server_local_canonical_without_flag_rejected',records=changed)
    add('server_local_threads_without_flag_rejected',env_overrides={'OMP_NUM_THREADS':'12'})
    fixture={'scope':'Pure independent policy fixtures; not a numerical or checkpoint-restore result',
        'source_config_text':text,'source_bundle_sha':SOURCE_BUNDLE,'source_binary_sha':manifest['binary_sha256'],
        'source_config_sha':SOURCE_HASH,'source_input_sha':manifest['input_manifest_sha256'],
        'target_input_sha':manifest['input_manifest_sha256'],'actual_locked_source_cpu':expected_cpu,
        'environment':env,'cases':cases}
    with args.output.open('x',encoding='utf-8') as stream:
        json.dump(fixture,stream,ensure_ascii=False)
    print(json.dumps({'fixture':str(args.output),'cases':len(cases),'positive':sum(c['expected'] for c in cases),
        'server_cases':sum(c['profile']=='server' for c in cases),'local_cases':sum(c['profile']=='local' for c in cases),
        'sha256':hashlib.sha256(args.output.read_bytes()).hexdigest(),'no_PEA_or_KF_executed':True}))


if __name__ == '__main__':
    main()
