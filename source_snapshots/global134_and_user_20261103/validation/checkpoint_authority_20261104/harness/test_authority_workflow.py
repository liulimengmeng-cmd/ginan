"""Pure conditional-workflow contracts; process creation is forbidden."""
import argparse
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('_authority_workflow_contracts',HERE/'workflow_authority_hour.py')
W=importlib.util.module_from_spec(spec);spec.loader.exec_module(W)

def plan():
    return dict(schema=W.SCHEMA,nodes=W.nodes(),execution_plans={'replay_e300':{},'formal':{}},
        original_source='/home/rx/GINAN/us20_global134_osb5h_20261031_v1/run/checkpoints/E29-20240717015930-e240-r32f626de9463',
        automatic_retry_authorised=False,resource_abort_thresholds_enabled=False,server_access_authorised=False,
        new_solver_started_by_prepare=False,formal_requires_actual_e300_restore=True)

def result(node):
    value=dict(schema='CHECKPOINT_AUTHORITY_AR_ACTION_RESULT_V1',passed=True,profile=node['profile'],
        action=node['action'],execution_plan_sha256='a'*64,original_source_objects_configuration_inputs_unchanged=True)
    count=60 if node['profile']=='replay_e300' else 120
    if node['action']=='validation':value.update(actual_fresh_PEA_restores=2,all_nine_section_audits=2,
        stages=[{},{}],finite_data_run_started=False)
    else:
        value.update(complete_epochs=count,complete_FLOAT_epochs=count,AR_enabled_epochs=count)
        if node['action']=='audit':
            captures=1 if node['profile']=='replay_e300' else 2
            value.update(atomic_product_commits=count,products=[{}]*count,
                regression_product_events=[{'passed':True}]*5,
                all_captured_checkpoints_actual_roundtrip_tested=True,actual_checkpoint_restores=[{}]*captures,
                checkpoints={'300':{}} if captures==1 else {'300':{},'360':{}})
    return value

class Workflow(unittest.TestCase):
    def test_fixed_six_nodes(self):
        self.assertEqual(W.nodes(),[dict(profile=p,action=a) for p in ('replay_e300','formal')
            for a in ('validation','formal','audit')]);W.validate_contract(plan())
    def test_no_reorder_skip_or_repeat(self):
        original=plan()
        for value in (original['nodes'][1:],original['nodes'][::-1],original['nodes']+[original['nodes'][0]]):
            bad=copy.deepcopy(original);bad['nodes']=value
            with self.assertRaises(RuntimeError):W.validate_contract(bad)
    def test_no_retry_resource_guard_or_server_scope(self):
        for field in ('automatic_retry_authorised','resource_abort_thresholds_enabled','server_access_authorised','new_solver_started_by_prepare'):
            bad=plan();bad[field]=True
            with self.assertRaises(RuntimeError):W.validate_contract(bad)
    def test_cannot_skip_actual_e300_restore(self):
        bad=plan();bad['formal_requires_actual_e300_restore']=False
        with self.assertRaises(RuntimeError):W.validate_contract(bad)
    def test_no_other_source_or_missing_plan(self):
        bad=plan();bad['original_source']+='/../other'
        with self.assertRaises(RuntimeError):W.validate_contract(bad)
        bad=plan();del bad['execution_plans']['replay_e300']
        with self.assertRaises(RuntimeError):W.validate_contract(bad)
    def test_exact_root_gate(self):
        gate=dict(schema=W.GATE_SCHEMA,approved=True,workflow_plan_sha256='a'*64,nodes=W.nodes(),
            automatic_retry_authorised=False,allow_formal_only_after_all_replay_gates_pass=True)
        W.validate_gate(gate,'a'*64)
        for key,value in (('approved',False),('workflow_plan_sha256','b'*64),('nodes',W.nodes()[1:]),
            ('automatic_retry_authorised',True),('allow_formal_only_after_all_replay_gates_pass',False)):
            bad=dict(gate);bad[key]=value
            with self.assertRaises(RuntimeError):W.validate_gate(bad,'a'*64)
    def test_predecessor_prefix_must_be_exact(self):
        done=[dict(node=node,passed=True) for node in W.nodes()[:3]];W.completed_prefix(done,3)
        for altered in (done[:2],done[::-1],done+[done[0]]):
            with self.assertRaises(RuntimeError):W.completed_prefix(altered,3)
        bad=copy.deepcopy(done);bad[2]['passed']=False
        with self.assertRaises(RuntimeError):W.completed_prefix(bad,3)
    def test_completed_receipts_remain_bound_to_original_pins(self):
        with tempfile.TemporaryDirectory() as d,patch.object(W,'kernel',return_value=None):
            root=Path(d);node=W.nodes()[0]
            who=dict(pid=10,birth_ticks=100,sid=10,uid=1000,executable='/usr/bin/python3.10')
            launch=dict(execution_plan_sha256='a'*64,action=node['action'],action_nonce='b'*32,controller=who)
            value=result(node);value.update(action_nonce='b'*32,controller_identity=who)
            launch_path=root/'launch.json';result_path=root/'result.json'
            launch_path.write_text(json.dumps(launch));result_path.write_text(json.dumps(value))
            completed=[dict(node=node,passed=True,result=W.identity(result_path),launch=W.identity(launch_path),
                exact_child_waitpid=dict(pid=10,birth_ticks=100,exit_code=0))]
            execution=dict(execution_plans={'replay_e300':{'sha256':'a'*64}})
            W.verify_completed_prefix(execution,completed,1)
            altered=dict(value,unreviewed_extra_field=True);result_path.write_text(json.dumps(altered))
            with self.assertRaises(RuntimeError):W.verify_completed_prefix(execution,completed,1)
            result_path.write_text(json.dumps(value));launch_path.write_text(json.dumps(dict(launch,action_nonce='c'*32)))
            with self.assertRaises(RuntimeError):W.verify_completed_prefix(execution,completed,1)
    def test_all_nodes_require_exact_result(self):
        for node in W.nodes():
            value=result(node);W.validate_result(value,node,'a'*64)
            for key,wrong in (('passed',False),('profile','wrong'),('action','wrong'),
                ('execution_plan_sha256','b'*64),('original_source_objects_configuration_inputs_unchanged',False)):
                bad=dict(value);bad[key]=wrong
                with self.assertRaises(RuntimeError):W.validate_result(bad,node,'a'*64)
    def test_envelope_is_not_actual_restore(self):
        for node in W.nodes():
            if node['action']!='validation':continue
            for field,value in (('actual_fresh_PEA_restores',0),('all_nine_section_audits',0),
                ('finite_data_run_started',True),('stages',[{}])):
                bad=result(node);bad[field]=value
                with self.assertRaises(RuntimeError):W.validate_result(bad,node,'a'*64)
    def test_all_actual_epochs_required(self):
        for node in W.nodes():
            if node['action']=='validation':continue
            for field in ('complete_epochs','complete_FLOAT_epochs','AR_enabled_epochs'):
                bad=result(node);bad[field]-=1
                with self.assertRaises(RuntimeError):W.validate_result(bad,node,'a'*64)
    def test_audit_requires_actual_checkpoint_restores_and_products(self):
        for node in W.nodes():
            if node['action']!='audit':continue
            for key,value in (('actual_checkpoint_restores',[]),('products',[]),('regression_product_events',[]),('checkpoints',{}),
                ('all_captured_checkpoints_actual_roundtrip_tested',False),('atomic_product_commits',0)):
                bad=result(node);bad[key]=value
                with self.assertRaises(RuntimeError):W.validate_result(bad,node,'a'*64)
    def test_pid_alone_never_establishes_owner(self):
        who=dict(pid=10,birth_ticks=100,sid=10,uid=1000,executable='/usr/bin/python3.10')
        self.assertTrue(W.same_owner(who,dict(who,state='R')))
        for key,value in (('birth_ticks',101),('sid',11),('uid',0),('executable','/other/python')):
            wrong=dict(who);wrong[key]=value;self.assertFalse(W.same_owner(who,wrong))
        self.assertFalse(W.same_owner(who,None))
    def test_duplicate_nonfinite_json_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            file=Path(d)/'a.json'
            for raw in ('{"a":1,"a":2}','{"a":NaN}','{"a":1e999}'):
                file.write_text(raw)
                with self.assertRaises(RuntimeError):W.load(file)
    def test_strict_pin_content(self):
        with tempfile.TemporaryDirectory() as d:
            file=Path(d)/'a';file.write_text('first');pin=W.identity(file);W.verify(pin)
            file.write_text('second')
            with self.assertRaises(RuntimeError):W.verify(pin)
    def test_exclusive_failure_does_not_claim_canonical_result(self):
        from types import SimpleNamespace
        fake=dict(pid=10,birth_ticks=100,sid=10,uid=1000,executable='/usr/bin/python3.10')
        with tempfile.TemporaryDirectory() as d,patch.object(W,'DIRECTORY',Path(d)),patch.object(W,'kernel',return_value=fake):
            args=SimpleNamespace(nonce='../invalid',plan_sha256='a'*64)
            with patch('sys.stderr'):
                value=W.exclusive_failure(args,RuntimeError('test'),'STARTUP')
            self.assertFalse((Path(d)/'RESULT.json').exists())
            self.assertFalse((Path(d)/'STARTED.json').exists())
            self.assertTrue(value['canonical_ownership_not_granted_by_this_receipt'])
            self.assertEqual(len(list(Path(d).glob('STARTUP--unverified-*.json'))),1)

def main():
    p=argparse.ArgumentParser();p.add_argument('--receipt',type=Path);args=p.parse_args()
    suite=unittest.defaultTestLoader.loadTestsFromModule(__import__(__name__))
    with patch('subprocess.Popen',side_effect=AssertionError('No real launch in pure tests')),\
         patch('subprocess.run',side_effect=AssertionError('No real launch in pure tests')):
        ran=unittest.TextTestRunner(verbosity=2).run(suite)
    receipt=dict(passed=ran.wasSuccessful(),tests_run=ran.testsRun,solver_processes_started=0,
        workflow=W.identity(HERE/'workflow_authority_hour.py'),test_source=W.identity(__file__),
        failures=len(ran.failures),errors=len(ran.errors),scientific_accuracy_accepted=False)
    if args.receipt:W.write_new(args.receipt,receipt)
    if not ran.wasSuccessful():raise SystemExit(1)

if __name__=='__main__':main()
