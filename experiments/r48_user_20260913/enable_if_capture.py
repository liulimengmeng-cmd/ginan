#!/usr/bin/env python3
"""Enable the factor buffer required by the current IF joint-noise callbacks."""
import hashlib,json
from pathlib import Path
import yaml
W=Path('/mnt/c/Users/rx/Documents/GINAN/r48_user_work_20260913_prefix')
s=json.loads((W/'experiment.json').read_text())
for c in s['cases']:
 cfg=yaml.safe_load(Path(c['config']).read_text());c['case']+='_capture';c['output']+='_capture'
 assert not Path(c['output']).exists()
 z=cfg['processing_options']['gnss_general']['zhang_pppar'];z['fixed_lag_factor_capture_shadow']=True;z['checkpoint_runtime_id']=c['case']
 cfg['outputs']['metadata']['config_description']=c['case']
 cfg['receiver_options']['global']['models']['ionospheric_components'].update(use_2nd_order=False,use_3rd_order=False)
 p=W/'config'/(c['case']+'.yaml');assert not p.exists();p.write_text(yaml.safe_dump(cfg,sort_keys=False))
 c['config']=str(p);c['config_sha256']=hashlib.sha256(p.read_bytes()).hexdigest()
s['configuration_correction']='IF joint-noise transition callbacks require the E18 factor capture buffer in the current frozen binary. Enabled equally in FLOAT and AR; old failed run preserved.'
p=W/'experiment_capture.json';assert not p.exists();p.write_text(json.dumps(s,indent=2))
print(p)
