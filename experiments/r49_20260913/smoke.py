from pathlib import Path
import subprocess,os,re,json,time,sys
a=Path(__file__).resolve().parent;spec=json.loads((a/'experiment.json').read_text())
suffix=sys.argv[1] if len(sys.argv)>1 else ''
assert re.fullmatch(r'[a-z0-9_]*',suffix)
smoke=spec['case']+'_checkpoint_smoke'+suffix
output=Path('/mnt/d/GINAN_R20/inputData/outputs')/smoke
assert not output.exists()
s=(a/'config/case.yaml').read_text().replace(spec['case'],smoke)
s=s.replace('end_epoch: 2024-07-17 00:30:00','end_epoch: 2024-07-17 00:00:00')
checkpoint='/mnt/d/GINAN_R20/r49_checkpoint_smoke_20260913'+suffix
s=s.replace(spec['checkpoints']['directory'],checkpoint)
s=re.sub(r'      checkpoint_capture_epochs:\n(?:        - .*\n)+','      checkpoint_capture_epochs:\n        - "2024-07-17 00:00:00"\n',s)
(a/'smoke.yaml').write_text(s)
files=['zhang_global_2024199_180_base.yaml','zhang_global_2024199_180_inputs.yaml',
 'zhang_global_2024199_180_product.yaml','zhang_global_2024199_180_e29_180network_product_1h.yaml',
 'zhang_global_2024199_180_e29_service_30s_1h.yaml']
cmd=[spec['source_root']+'/bin/pea','-q','-y',*[str(a/'config'/f) for f in files],str(a/'smoke.yaml'),'-d',smoke]
start=time.time()
with (a/'smoke.log').open('w') as f:
 rc=subprocess.run(cmd,cwd='/mnt/d/GINAN_R20/inputData',env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'4','MKL_NUM_THREADS':'1','ZHANG_R49_FUSION_WEIGHT':'0.10'},stdout=f,stderr=subprocess.STDOUT).returncode
out={'status':rc,'seconds':time.time()-start,'output':str(output),'checkpoint_directory':checkpoint}
(a/'smoke_result.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out));assert rc==0
