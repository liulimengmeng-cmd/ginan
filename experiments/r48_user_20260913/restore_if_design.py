#!/usr/bin/env python3
"""Retain the STEC design columns required to form the IF observation combination."""
import hashlib,json
from pathlib import Path
import yaml
W=Path('/mnt/c/Users/rx/Documents/GINAN/r48_user_work_20260913_prefix')
s=json.loads((W/'experiment_capture.json').read_text())
for c in s['cases']:
 cfg=yaml.safe_load(Path(c['config']).read_text());c['case']+='_if';c['output']+='_if'
 assert not Path(c['output']).exists()
 cfg['estimation_parameters']['receivers']['global']['ion_stec']['estimated']=[True]
 cfg['processing_options']['gnss_general']['zhang_pppar']['checkpoint_runtime_id']=c['case']
 cfg['outputs']['metadata']['config_description']=c['case']
 p=W/'config'/(c['case']+'.yaml');assert not p.exists();p.write_text(yaml.safe_dump(cfg,sort_keys=False))
 c['config']=str(p);c['config_sha256']=hashlib.sha256(p.read_bytes()).hexdigest()
s['if_design_correction']='Restore ion_stec.estimated=true required by Ginan to construct IF combinations from ionosphere design columns. Both controls identical. Earlier runs are failed preflight cases.'
p=W/'experiment_if.json';assert not p.exists();p.write_text(json.dumps(s,indent=2))
print(p)
