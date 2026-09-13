#!/usr/bin/env python3
"""Prepare UDUC paired controls from the exact completed IF experiment."""
import hashlib,json,subprocess
from pathlib import Path
import yaml
HERE=Path(__file__).resolve().parent
W=Path('/mnt/c/Users/rx/Documents/GINAN/r48_uduc_work_20260913')
OLD=Path('/mnt/c/Users/rx/Documents/GINAN/r48_user_work_20260913_prefix')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 W.mkdir(exist_ok=False);(W/'config').mkdir()
 s=json.loads((OLD/'experiment_if.json').read_text())
 assert sha(s['binary'])==s['binary_sha256']
 for v in s['snapshot'].values():assert sha(v['path'])==v['sha256']
 for st in s['stations']:assert sha(st['path'])==st['sha256'] and not st['in_server_network']
 for c in s['cases']:
  assert sha(c['config'])==c['config_sha256']
  cfg=yaml.safe_load(Path(c['config']).read_text())
  c['if_case']=c['case'];c['if_config']=c['config'];c['if_output']=c['output']
  c['case']=f"r48_uduc_20260913_{c['station']}_{c['mode']}"
  c['output']=str(Path(s['working_directory'])/'outputs'/c['case']);assert not Path(c['output']).exists()
  z=cfg['processing_options']['gnss_general']['zhang_pppar'];z['integer_strategy']='CANONICAL_USER_SD_WL_L1';z['checkpoint_runtime_id']=c['case']
  assert z['fixed_lag_factor_capture_shadow'] and z['user_use_full_product_covariance']
  io=cfg['processing_options']['ppp_filter']['ionospheric_components'];io.update(common_ionosphere=True,use_if_combo=False,use_gf_combo=False,corr_mode='estimate')
  assert cfg['estimation_parameters']['receivers']['global']['ion_stec']['estimated']==[True]
  cfg['outputs']['metadata']['config_description']=c['case']
  p=W/'config'/(c['case']+'.yaml');p.write_text(yaml.safe_dump(cfg,sort_keys=False));c['config']=str(p);c['config_sha256']=sha(p)
 s.update(experiment='R48 UDUC user FLOAT versus AR paired control',work=str(W),if_comparison_work=str(OLD),experiment_commit=subprocess.check_output(['git','-C',str(HERE.parent.parent),'rev-parse','HEAD'],text=True).strip(),claims='Same 44 completed PRODUCT_FIXED epochs and five held-out users as IF. Full product covariance and factor replay audit retained. UDUC explicitly estimates STEC without extra ionosphere constraints. WL-only fixing must not be labelled full dual-frequency AR. Upstream DGEMV remains unresolved.')
 (W/'experiment.json').write_text(json.dumps(s,indent=2));print(json.dumps({'work':str(W),'cases':len(s['cases']),'epochs':s['expected_epochs']},indent=2))
if __name__=='__main__':main()
