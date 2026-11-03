"""Pure finite-plan/lifecycle contracts. Every subprocess launch is forbidden."""
import argparse
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
def module(name,file):
    spec=importlib.util.spec_from_file_location(name,HERE/file);value=importlib.util.module_from_spec(spec);spec.loader.exec_module(value);return value
P=module('_datum_plan_test','prepare_authority_hour.py')
R=module('_datum_run_test','run_authority_hour.py')
RAW=P.host_path(P.ART.parent/'ar_checkpoint_codec_fix_20261103_01/checkpoint_inputs/source_e240_manifest.json').read_bytes()
E=R.configure('formal');C=E.C

def plan(profile='formal'):return P.build_plan(RAW,'a'*64,P.AUDITOR,profile)
def source():return P.P.decode(json.loads(RAW)['config_text'].encode())
def complete():
    marker=C.Markers(HERE/'NO_RUN_CREATED',HERE/'NO_RUNNER_CREATED',False);label='/fixture/Network.TRACE'
    marker.ingest(label,'R51_VALIDATION_POLICY ratio_only=1 ratio_threshold=2 nis_gate=0 integer_product_nis_evaluation=SKIPPED_RATIO_ONLY')
    marker.ingest(label,'US20_E240_NATIVE_DECODED_RESTORE_ADMITTED')
    for i,t in enumerate(R.epoch_grid('formal')):
        n=1000*i+100
        marker.ingest(label,f'ZHANG_FLOAT_PHASE_BEGIN phase=FLOAT_TOTAL epoch="{t}" monotonic_ns={n}')
        for boundary,status in [('B13','COMMITTED'),('B16','AUTHORITATIVE_COMMIT')]:
            marker.ingest(label,f'ZHANG_KF_MEASUREMENT_TRANSACTION time={t} boundary={boundary} status={status} transaction_id={t}|/PPP|1 model_generation=3 solve_generation=4 final_commit_sequence=5 qc_converged=1')
        marker.ingest(label,f'ZHANG_FLOAT_PHASE_TIMER phase=FLOAT_TOTAL epoch="{t}" monotonic_ns={n+100}')
        marker.ingest(label,f'ZHANG_P0_RESOURCE phase=AR_CALL epoch="{t}" boundary=BEGIN monotonic_ns={n+110}')
        marker.ingest(label,f'R51_AR_SCHEDULE time={t} phase=AR_ENABLED float_state_reset=0')
        marker.ingest(label,f'ZHANG_R51_EPOCH_ONLY_AR time={t} epoch_transaction_started=1 prior_integer_reuse=0 float_filter_preserved=1 physical_arcs_preserved=1')
        marker.ingest(label,f'ZHANG_PER_EPOCH_AR_OWNER time={t} authoritative_selected=1 transactional=1 legacy_hold=0 authoritative_object=0x101 selected_object=0x101 augmented_object=0x202')
        marker.ingest(label,f'ZHANG_P0_RESOURCE phase=AR_CALL epoch="{t}" boundary=END monotonic_ns={n+210}')
    return marker
def sections(ppp=6,gec=1):
    return {'sections':[dict(name=n,schema_version=ppp if n=='zhang_ppp_ar_runtime' else gec if n=='GEC_AR_RUNTIME_V1' else v,
        payload_bytes=20,sha256='a'*64) for n,v in C.SECTIONS.items()]}

class Plans(unittest.TestCase):
    def test_hour_exact120(self):
        p,c,o=plan();self.assertEqual((p['target']['first_epoch'],p['target']['last_epoch'],p['target']['new_epochs']),(241,360,120))
        self.assertEqual(p['target']['end_gpst'],'2024-07-17 02:59:30');self.assertEqual(len(R.epoch_grid('formal')),120)
    def test_no_warmup(self):
        p,_,_=plan();self.assertEqual(p['acceptance']['expected_FLOAT_only_epochs'],0)
        self.assertEqual(p['acceptance']['required_actual_AR_epochs'],120)
        self.assertTrue(all(R.phase_for('formal',t)=='AR' for t in R.epoch_grid('formal')))
    def test_captures(self):
        p,_,_=plan();self.assertEqual(p['target']['captures'],[dict(epoch=300,gpst='2024-07-17 02:29:30'),dict(epoch=360,gpst='2024-07-17 02:59:30')])
    def test_probe_is_independent(self):
        p,_,_=plan('replay_e300');self.assertEqual(p['target']['new_epochs'],60);self.assertEqual(p['target']['last_epoch'],300)
        self.assertEqual(p['target']['run_root'],P.ROOT+'/validation_e300/run')
        self.assertEqual([r['epoch'] for r in p['target']['captures']],[300])
    def test_source_pin(self):
        p,_,_=plan();self.assertEqual((p['source']['completed_epoch'],p['source']['state_dimension']),(240,7274))
        self.assertEqual(p['source']['bundle_sha256'],P.SOURCE_BUNDLE_SHA);self.assertEqual(len(p['source']['sections']),9)
    def test149_argv_preserved(self):
        old=source();argv=[old[f'command_argument_{i}'] for i in range(1,150)]
        for profile in P.PROFILES:
            p,c,_=plan(profile);actual=p['real_run_after_separate_root_gate']['argv'][1:]
            self.assertEqual(actual[:10]+actual[11:],argv);self.assertEqual(len(actual),150)
            self.assertEqual(P.P.decode(c)['config_file_count'],'11')
    def test_all10_config_bytes_preserved(self):
        old=source();_,c,_=plan();new=P.P.decode(c)
        for i in range(10):
            for suffix in ('path','bytes','sha256','size'):self.assertEqual(new[f'config_file_{i}_{suffix}'],old[f'config_file_{i}_{suffix}'])
    def test_exact_overlay(self):
        _,_,o=plan();self.assertIn(b"end_epoch: '2024-07-17 02:59:30'\n    max_epochs: 360",o)
        self.assertIn(b"config_description: 'Checkpoint authority repair direct AR1h'",o)
        self.assertIn(P.ROOT.encode(),o)
    def test_four_source_gates(self):
        p,_,_=plan();g=p['ordered_pre_run_gates'];self.assertEqual(len(g),4)
        self.assertEqual(g[0]['environment']['ZHANG_GEC134_RESTORE_ONLY'],'1')
        self.assertEqual(g[0]['environment']['ZHANG_E29_RESTORE_BUNDLE'],P.SOURCE)
        self.assertNotEqual(g[2]['environment']['ZHANG_E29_RESTORE_BUNDLE'],P.SOURCE)
        self.assertEqual([g[i]['argv'][-1] for i in (1,3)],['source-migration','same-producer'])
    def test_two_own_real_restores(self):
        p,_,_=plan();g=p['ordered_post_run_gates'];self.assertEqual(len(g),4)
        self.assertEqual([g[i]['completed_epoch'] for i in (0,2)],[300,360])
        self.assertTrue(all(g[i]['environment']['ZHANG_GEC134_RESTORE_ONLY']=='1' for i in (0,2)))
    def test_policy_unchanged(self):
        p,_,_=plan();a=p['acceptance'];self.assertEqual((a['ratio'],a['joint_NIS'],a['epoch_independent_AR'],a['WSL_GiB']),(2,False,True,20))
        self.assertEqual(p['resource_policy'],dict(new_resource_abort_thresholds=False,observational_samples_only=True))
        env=p['real_run_after_separate_root_gate']['environment']
        for k,v in {'OMP_NUM_THREADS':'12','OPENBLAS_NUM_THREADS':'1','ZHANG_FLOAT_MATRIX_THREADS':'6',
            'ZHANG_R51_RATIO_THRESHOLD':'2','ZHANG_R51_EPOCH_ONLY_AR':'1','ZHANG_REGIONAL_JOINT_EXPORT':'0',
            'ZHANG_REGIONAL_TROPOSPHERE':'OFF'}.items():self.assertEqual(env[k],v)
    def test_no_old_product_equality_or_probe_gate(self):
        p,_,_=plan();self.assertNotIn('matched_original_e241_product_gate',p)
        self.assertNotIn('reversible_first_pass_archive',p);self.assertFalse(p['acceptance']['old_algorithm_product_value_equality_required'])
    def test_bad_authority_rejected(self):
        for raw,target,auditor,profile in ((RAW,'a'*64,P.AUDITOR,'bad'),(b'{}','a'*64,P.AUDITOR,'formal'),
            (RAW,P.SOURCE_BINARY_SHA,P.AUDITOR,'formal'),(RAW,'A'*64,P.AUDITOR,'formal'),(RAW,'a'*64,'/other/audit','formal')):
            with self.assertRaises((ValueError,RuntimeError)):P.build_plan(raw,target,auditor,profile)

class Contracts(unittest.TestCase):
    def test_allowed_source_formats(self):
        for ppp in (6,7):
            for gec in (1,2):self.assertEqual(len(R.section_registry(sections(ppp,gec))),9)
    def test_bad_versions(self):
        for ppp,gec in ((5,1),(8,1),(6,0),(6,3),(True,1),(6,True)):
            with self.assertRaises(RuntimeError):R.section_registry(sections(ppp,gec))
    def test_missing_duplicate_sections(self):
        for mutate in (lambda x:x['sections'].pop(),lambda x:x['sections'].append(x['sections'][0])):
            s=sections();mutate(s)
            with self.assertRaises(RuntimeError):R.section_registry(s)
    def test_nonpositive_payload(self):
        s=sections();s['sections'][0]['payload_bytes']=0
        with self.assertRaises(RuntimeError):R.section_registry(s)
    def test_all120_real_markers(self):
        m=complete();rows=R.audit_epochs(m,'formal');self.assertEqual(len(rows),120)
        self.assertTrue(all(r['AR_enabled_epoch'] for r in rows));self.assertEqual(m.summary()['AR_enabled_epochs'],120)
        self.assertEqual(m.summary()['epochs_with_R51_SEARCH_ENTERED'],0)
    def test_enabled_is_not_search_count(self):
        m=complete();m.ingest('/fixture/Network.TRACE',f'R51_SEARCH_ENTERED time={R.epoch_grid("formal")[0]}')
        self.assertEqual(m.summary()['epochs_with_R51_SEARCH_ENTERED'],1);self.assertEqual(len(R.audit_epochs(m,'formal')),120)
    def test_missing_epoch(self):
        m=complete();m.epochs.pop(R.epoch_grid('formal')[0])
        with self.assertRaises(RuntimeError):R.audit_epochs(m,'formal')
    def test_incomplete_float(self):
        m=complete();m.epochs[R.epoch_grid('formal')[0]].pop('float_end_ns')
        with self.assertRaises(RuntimeError):R.audit_epochs(m,'formal')
    def test_authoritative_float_rollback_rejected(self):
        t=R.epoch_grid('formal')[0];m=complete()
        with self.assertRaises(RuntimeError):
            m.ingest('/fixture/Network.TRACE',f'ZHANG_KF_MEASUREMENT_TRANSACTION time={t} boundary=B13 status=ROLLED_BACK transaction_id={t}|/PPP|1 qc_converged=0')
            R.audit_epochs(m,'formal')
    def test_plain_ar_rollback_not_float_failure(self):
        m=complete();m.ingest('/fixture/Network.TRACE','ZHANG_FIXED_SUBTRANSACTION stage=L1 status=ROLLED_BACK')
        self.assertEqual(len(R.audit_epochs(m,'formal')),120)
    def test_warmup_rejected(self):
        m=complete()
        with self.assertRaises(RuntimeError):m.ingest('/fixture/Network.TRACE',f'R51_AR_SCHEDULE time={R.epoch_grid("formal")[0]} phase=FLOAT_WARMUP search_calls=0 history_applied=0 authoritative_feedback=0')
    def test_unknown_epoch_rejected(self):
        with self.assertRaises(RuntimeError):R.phase_for('formal','2024-07-17 03:00:00')
    def test_gate_is_exact_and_single_use_policy(self):
        gate=dict(schema=R.GATE_SCHEMA,action='validation',profile='replay_e300',approved=True,execution_plan_sha256='a'*64,automatic_retry_authorised=False)
        with patch.object(R,'PROFILE','replay_e300'):
            R.validate_gate('validation',gate,'a'*64)
            for k,v in (('approved',False),('profile','formal'),('automatic_retry_authorised',True),('execution_plan_sha256','b'*64)):
                bad=dict(gate);bad[k]=v
                with self.assertRaises(RuntimeError):R.validate_gate('validation',bad,'a'*64)
    def test_formal_requires_real_replay(self):
        gate=dict(schema=R.GATE_SCHEMA,action='validation',profile='formal',approved=True,execution_plan_sha256='a'*64,automatic_retry_authorised=False)
        with patch.object(R,'PROFILE','formal'):
            with self.assertRaises(RuntimeError):R.validate_gate('validation',gate,'a'*64)
    def test_replay_acceptance_is_candidate_bound_and_complete(self):
        with tempfile.TemporaryDirectory() as d,patch.object(R,'ART',Path(d)):
            root=Path(d);replay=root/'controller_replay_e300_01';formal=root/'controller_formal_01'
            replay.mkdir();formal.mkdir()
            plan=dict(qualification={'sha256':'a'*64},science=dict(source={'bundle':'s'},target={'binary_sha256':'b'*64}))
            for target in (replay,formal):(target/'EXECUTION_PLAN.json').write_text(json.dumps(plan))
            result=dict(passed=True,profile='replay_e300',execution_plan_sha256=R.sha(replay/'EXECUTION_PLAN.json'),
                complete_epochs=60,atomic_product_commits=60,complete_FLOAT_epochs=60,AR_enabled_epochs=60,
                all_captured_checkpoints_actual_roundtrip_tested=True,actual_checkpoint_restores=[{}],
                products=[{}]*60,regression_product_events=[{'passed':True}]*5,
                checkpoints={'300':{}},original_source_objects_configuration_inputs_unchanged=True)
            target=replay/'audit_RESULT.json';target.write_text(json.dumps(result));R.replay_acceptance(R.identity(target))
            for key,value in (('passed',False),('complete_epochs',59),('actual_checkpoint_restores',[]),
                ('products',[{}]*59),('all_captured_checkpoints_actual_roundtrip_tested',False),('checkpoints',{})):
                bad=dict(result);bad[key]=value;target.write_text(json.dumps(bad))
                with self.assertRaises(RuntimeError):R.replay_acceptance(R.identity(target))
            target.write_text(json.dumps(result));wrong=copy.deepcopy(plan);wrong['science']['target']['binary_sha256']='c'*64
            (formal/'EXECUTION_PLAN.json').write_text(json.dumps(wrong))
            with self.assertRaises(RuntimeError):R.replay_acceptance(R.identity(target))
    def test_pending_gps_does_not_cancel_other_systems(self):
        commit=dict(epoch='2024-07-17 02:27:30',route='ENHANCED',ar=True,
            systems={system:dict(pppar_usable_phase_rows=count) for system,count in [('GPS',0),('GAL',4),('BDS',4)]})
        row=dict(satellite='G03',ar_valid='0',pppar_usable='0',integer_valid='0',integer_transport_valid='0',dual_frequency_ar_valid='0')
        value=R.validate_event_product(commit,[row]);self.assertEqual(value['dual_usable'],dict(GPS=0,GAL=2,BDS=2))
        for key in ('ar_valid','pppar_usable','integer_valid','integer_transport_valid','dual_frequency_ar_valid'):
            bad=dict(row);bad[key]='1'
            with self.assertRaises(RuntimeError):R.validate_event_product(commit,[bad])
        for system in ('GAL','BDS'):
            bad=copy.deepcopy(commit);bad['systems'][system]['pppar_usable_phase_rows']=2
            with self.assertRaises(RuntimeError):R.validate_event_product(bad,[row])
        bad=copy.deepcopy(commit);bad['systems']['GPS']['pppar_usable_phase_rows']=2
        with self.assertRaises(RuntimeError):R.validate_event_product(bad,[row])
        bad=dict(commit,route='FLOAT_ONLY',ar=False)
        with self.assertRaises(RuntimeError):R.validate_event_product(bad,[row])
    def test_failed_run_metadata_does_not_follow_input_links(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)/'failed';root.mkdir();outside=Path(d)/'inputs';outside.mkdir()
            (root/'small').write_text('first');(outside/'raw').write_text('observation')
            (root/'input_link').symlink_to(outside,target_is_directory=True)
            value=R.failed_run_metadata([root]);entries=value[str(root)]
            self.assertEqual(entries['input_link']['kind'],'link')
            self.assertNotIn('input_link/raw',entries)
            (root/'small').write_text('changed')
            self.assertNotEqual(value,R.failed_run_metadata([root]))
    def test_exact_source_schema_report(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'checkpoint_manifest.json').write_text(json.dumps(sections(6)))
            report=dict(passed=True,audit_kind='same-producer',failed_fields=[],estimation_executed_by_auditor=False,
                whole_process_restore_must_be_proven_by_PEAs_native_receipt=True,scientific_accuracy_accepted=False,
                checks={f'core.{i}':True for i in range(33)},state_dimension=2,P_entries_raw_IEEE=4,
                sections=[dict(name=x['name'],schema=x['schema_version'],entire_persistent_payload_equal=True,
                    declared_ready_workset_boundary_policy_only=x['name']=='pea.receiver_runtime.v1') for x in sections(6)['sections']])
            path=root/'audit.json';path.write_text(json.dumps(report));R.validate_nine_report(path,'same-producer',root)
            next(x for x in report['sections'] if x['name']=='zhang_ppp_ar_runtime')['schema']=7
            path.write_text(json.dumps(report))
            with self.assertRaises(RuntimeError):R.validate_nine_report(path,'same-producer',root)
    def test_strict_json(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'a.json'
            for text in ('{"a":1,"a":2}','{"a":NaN}','{"a":1e999}'):
                p.write_text(text)
                with self.assertRaises(RuntimeError):R.strict_load(p)
    def test_adapted_run_has_no_old_epoch_counts(self):
        def values(code):
            for x in code.co_consts:
                if hasattr(x,'co_consts'):yield from values(x)
                else:yield x
        self.assertFalse(any(type(x) is int and x in (600,840) for x in values(E.run.__code__)))

def main():
    p=argparse.ArgumentParser();p.add_argument('--receipt',type=Path);args=p.parse_args()
    suite=unittest.defaultTestLoader.loadTestsFromModule(__import__(__name__))
    with (patch('subprocess.Popen',side_effect=AssertionError('No process launch in pure tests')),
          patch('subprocess.run',side_effect=AssertionError('No process launch in pure tests'))):
        result=unittest.TextTestRunner(verbosity=2).run(suite)
    receipt=dict(passed=result.wasSuccessful(),tests_run=result.testsRun,solver_processes_started=0,
        controller=R.identity(HERE/'run_authority_hour.py'),preparer=R.identity(HERE/'prepare_authority_hour.py'),test_source=R.identity(__file__),
        failed=len(result.failures),errors=len(result.errors),scientific_accuracy_accepted=False)
    if args.receipt:R.write_new(args.receipt,receipt)
    if not result.wasSuccessful():raise SystemExit(1)

if __name__=='__main__':main()
