#!/usr/bin/env python3
"""Freeze the R48 completed prefix and prepare five independent paired users."""
import csv, hashlib, json, math, re, shutil, subprocess
from collections import Counter, defaultdict
from pathlib import Path
import yaml
import numpy as np
HERE=Path(__file__).resolve().parent
WORK=Path('/mnt/c/Users/rx/Documents/GINAN/r48_user_work_20260913')
DATA=Path('/mnt/d/GINAN_R20/inputData')
SOURCE=HERE.parent.parent
SERVER=Path('/mnt/c/Users/rx/Documents/GINAN/r48fix1_work_20260913')
BINARY=Path('/home/rx/GINAN/frozen-r48fix1-active-basis-20260913/bin/pea')
FIRST=1405209600
LAST=FIRST+1470
END='2024-07-17 00:24:30'

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
 return h.hexdigest()

def merge(a,b):
 for k,v in b.items():
  if isinstance(v,dict) and isinstance(a.get(k),dict):merge(a[k],v)
  else:a[k]=v
 return a

def main():
 WORK.mkdir(exist_ok=False)
 (WORK/'products').mkdir();(WORK/'config').mkdir();(WORK/'evidence').mkdir()
 server=json.loads((SERVER/'experiment.json').read_text())
 raw=DATA/'outputs'/server['case']
 ck=Path(server['checkpoints']['directory'])/'E29-20240717002430-e50-r63ba52d74ee8/checkpoint_manifest.json'
 cm=json.loads(ck.read_text());assert cm['epoch']==END and (ck.parent/'checkpoint.bundle').is_file()
 assert sha(BINARY)==cm['binary_sha256']
 shutil.copy2(ck,WORK/'evidence/server_checkpoint_manifest.json')
 # The product and covariance prefixes were finalized before this checkpoint.
 snapshot={}
 for name in ['zhang_internal_products.csv','zhang_internal_product_covariance.csv']:
  src=raw/name;dst=WORK/'products'/name
  counts=Counter()
  with src.open() as fi,dst.open('w',newline='') as fo:
   rd=csv.DictReader(fi);wr=csv.DictWriter(fo,fieldnames=rd.fieldnames);wr.writeheader()
   for row in rd:
    epoch=int(float(row['gpst_seconds']))
    if epoch>LAST:continue
    if row['solution']!='PRODUCT_FIXED':continue
    assert FIRST<=epoch<=LAST and (epoch-FIRST)%30==0
    assert None not in row and all(v is not None for v in row.values())
    wr.writerow(row);counts[epoch]+=1
  assert set(counts)==set(range(FIRST,LAST+1,30))
  snapshot[name]={'source':str(src),'path':str(dst),'sha256':sha(dst),'rows_per_epoch':dict(sorted(counts.items()))}
 # Audit complete upper-triangle covariance and finite formal product fields.
 grouped=defaultdict(list)
 with (WORK/'products/zhang_internal_product_covariance.csv').open() as f:
  for r in csv.DictReader(f):grouped[int(r['gpst_seconds'])].append(r)
 cov_audit=[]
 for epoch,rows in sorted(grouped.items()):
  key=lambda r,s:(r[s+'_satellite'],r[s+'_parameter'],r[s+'_observable'])
  keys=sorted({key(r,s) for r in rows for s in ['row','column']});idx={k:i for i,k in enumerate(keys)};n=len(keys)
  seen=set();matrix=np.zeros((n,n))
  for r in rows:
   i,j=idx[key(r,'row')],idx[key(r,'column')];pair=tuple(sorted((i,j)))
   assert pair not in seen;seen.add(pair)
   v=float(r['covariance_m2']);assert math.isfinite(v);matrix[i,j]=matrix[j,i]=v
  assert len(seen)==n*(n+1)//2
  eig=np.linalg.eigvalsh(matrix);assert eig.min()>=-1e-8*max(1.,eig.max())
  cov_audit.append({'epoch':epoch,'dimension':n,'rows':len(rows),'min_eigenvalue':float(eig.min()),'max_eigenvalue':float(eig.max())})
 (WORK/'evidence/covariance_audit.json').write_text(json.dumps(cov_audit,indent=2))
 with (HERE/'validation_manifest.csv').open() as f: stations=list(csv.DictReader(f))[:5]
 net=yaml.safe_load((SERVER/'config/zhang_global_2024199_180_inputs.yaml').read_text())['inputs']['gnss_observations']['rnx_inputs']
 assert len(net)==180
 for station in stations:
  assert station['filename'] not in net and all(not x.startswith(station['station']) for x in net)
  obs=DATA/'data'/station['filename'];assert obs.is_file()
  station.update(path=str(obs),sha256=sha(obs),in_server_network=False)
 snx=DATA/'products/IGS0OPSSNX_20241960000_07D_07D_CRD.SNX'
 st=snx.read_text(); reference=[]
 for station in stations:
  ss=station['station'];coords={}
  for line in st.splitlines():
   fields=line.split()
   if len(fields)>8 and fields[1] in ['STAX','STAY','STAZ'] and fields[2]==ss:coords[fields[1]]=float(fields[8])
  assert len(coords)==3,ss
  reference.append({'station':ss,'coordinates_m':coords})
  for fn in ['OLOAD_ZHANG_2024199_180.BLQ','ALOAD_ZHANG_2024199_180.BLQ']:
   assert re.search(r'^\s*'+ss+r'\s*$',(DATA/'products/tables'/fn).read_text(),re.M)
 config={}
 for file in ['base.yaml','user_base.yaml','user_if.yaml']:merge(config,yaml.safe_load((HERE/file).read_text()))
 zp=config['processing_options']['gnss_general']['zhang_pppar']
 zp.update(product_filename=str(WORK/'products/zhang_internal_products.csv'),product_covariance_filename=str(WORK/'products/zhang_internal_product_covariance.csv'),product_solution='PRODUCT_FIXED',user_adapter=True,output_products=False,deterministic_checkpoint=False,checkpoint_runtime_id='R48-USER',canonical_user_target_feedback=True)
 config['processing_options']['epoch_control'].update(start_epoch='2024-07-17 00:00:00',end_epoch=END,epoch_interval=30,wait_next_epoch=3600)
 # IF processing does not estimate an unused STEC state.
 config['estimation_parameters']['receivers']['global']['ion_stec']['estimated']=[False]
 config['receiver_options']['global']['exclude']=False
 config['outputs']['outputs_root']='./outputs/<CONFIG>'
 cases=[]
 for station in stations:
  for mode in ['float','ar']:
   cfg=json.loads(json.dumps(config,default=str))
   case=f"r48_user_20260913_{station['station']}_{mode}"
   cfg['outputs']['metadata']['config_description']=case
   cfg['processing_options']['gnss_general']['zhang_pppar']['checkpoint_runtime_id']=case
   if mode=='float':
    cfg['processing_options']['ambiguity_resolution']['mode']='OFF'
    cfg['processing_options']['gnss_general']['zhang_pppar']['canonical_user_target_feedback']=False
   path=WORK/'config'/f'{case}.yaml';path.write_text(yaml.safe_dump(cfg,sort_keys=False))
   assert not (DATA/'outputs'/case).exists()
   cases.append({'station':station['station'],'mode':mode,'case':case,'config':str(path),'config_sha256':sha(path),'observation':station['path'],'output':str(DATA/'outputs'/case)})
 spec={'experiment':'R48 independent user PPP-AR paired prefix','work':str(WORK),'source_root':str(SOURCE),'server_source_commit':'94af7df19c5e5a153b5d8fb799c78a919e28892f','experiment_commit':subprocess.check_output(['git','-C',str(SOURCE),'rev-parse','HEAD'],text=True).strip(),'binary':str(BINARY),'binary_sha256':sha(BINARY),'working_directory':str(DATA),'first_epoch':FIRST,'last_epoch':LAST,'end':END,'expected_epochs':50,'stations':stations,'reference_sinex':{'path':str(snx),'sha256':sha(snx),'reference':reference},'snapshot':snapshot,'cases':cases,'threads':{'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'},'claims':'Independent of the 180-station product estimation; ENU is relative to IGS weekly SINEX apriori, not external ambiguity truth. Upstream DGEMV count 5 remains unresolved. Prefix only, not a full 30-minute or long-term validation.'}
 (WORK/'experiment.json').write_text(json.dumps(spec,indent=2))
 for f in (WORK/'products').iterdir():f.chmod(0o444)
 print(json.dumps({'work':str(WORK),'cases':len(cases),'epochs':50,'end':END,'covariance_epochs':len(cov_audit)},indent=2))
if __name__=='__main__':main()
