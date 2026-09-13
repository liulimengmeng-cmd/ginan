from pathlib import Path
import yaml,json,hashlib,datetime
w=Path('/mnt/c/Users/rx/Documents/GINAN/r48_user_model_work_20260914');s=json.loads((w/'experiment.json').read_text());cases=[]
for model in ['IF','UDUC']:
 for station in ['DYNG','NICO','BREW','MARS','JPLM']:
  base=Path('/mnt/c/Users/rx/Documents/GINAN/r48_uduc_work_20260913/config')/f'r48_uduc_20260913_{station}_ar.yaml' if model=='UDUC' else Path('/mnt/c/Users/rx/Documents/GINAN/r48_user_work_20260913_prefix/config')/f'r48_user_20260913_{station}_ar_capture_if.yaml'
  original=yaml.safe_load(base.read_text())
  for mode in ['float','round_shadow','round_feedback','ils_shadow','ils_feedback']:
   c=yaml.safe_load(base.read_text());name=f'r48u_20260914_{station}_{model}_{mode}'
   c['outputs']['metadata']['config_description']=name;c['outputs']['trace']['level']=4
   c['processing_options']['epoch_control']['end_epoch']='2024-07-17 00:30:00'
   c['processing_options']['ppp_filter']['joseph_stabilisation']=True
   z=c['processing_options']['gnss_general']['zhang_pppar'];z['product_filename']=str(w/'products/zhang_internal_products.csv');z['product_covariance_filename']=str(w/'products/zhang_internal_product_covariance.csv');z['checkpoint_runtime_id']=name;z['canonical_user_target_feedback']=mode.endswith('feedback')
   c['processing_options']['ambiguity_resolution']['mode']='OFF' if mode=='float' else 'LAMBDA'
   p=w/'config'/f'{name}.yaml';assert not p.exists();p.write_text(yaml.safe_dump(c,sort_keys=False))
   assert z['canonical_user_target_max_perr']==.001 and z['held_constraint_nis_alpha']==1e-6 and z['user_use_full_product_covariance']
   assert c['processing_options']['ambiguity_resolution']['solution_ratio_threshold']==3
   cases.append({'case':name,'station':station,'model':model,'mode':mode,'config':str(p),'output':'/mnt/d/GINAN_R20/inputData/outputs/'+name,'observation':next(x['path'] for x in s['stations'] if x['station']==station),'config_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'joint_ils':mode.startswith('ils')})
s['cases']=cases;s['joseph_stabilisation']=True;s['named_joint_mode']='LAMBDA_ALT';s['general_mixed_rows']='logged; only exact recovered named targets may pass original feedback gates';s['remaining_model_limits']=['single selected component reference remains in existing user implementation','MW physical and product eligibility lifecycle remains baseline; diagnostics retained','cross-epoch product covariance unavailable','no independent integer truth'];(w/'experiment.json').write_text(json.dumps(s,indent=2)+'\n')
record={'user_request':'Run R48 user AR tests concurrently with R49','started':'2026-09-14T00:11:42+08:00','user_threads':1,'priority':'nice 15; builds nice 19 serial','pure_r49_speed_comparison_valid_during_overlap':False,'server_binary_and_config_changed':False}
Path('/mnt/c/Users/rx/Documents/GINAN/r49_work_20260913/concurrent_user_tests_20260914.json').write_text(json.dumps(record,indent=2)+'\n')
print(len(cases),'matched full-window cases prepared')
