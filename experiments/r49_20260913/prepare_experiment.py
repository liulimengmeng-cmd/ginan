from pathlib import Path
import json,shutil,hashlib,datetime
w=Path(__file__).resolve().parent
old=w.parent/'r48fix1_work_20260913'
spec=json.loads((old/'experiment.json').read_text())
previous=spec['case'];case='zhang_p3_exp3_G3_r49_condition_domain_2024199_180_003000_20260913'
source='/home/rx/GINAN/r49-condition-domain-20260913'
frozen='/home/rx/GINAN/frozen-r49-condition-domain-20260913'
checkpoint='/mnt/d/GINAN_R20/r49_checkpoints_20260913'
spec['case']=case;spec['source_root']=source;spec['frozen_root']=frozen
spec['baseline_source_root']='/home/rx/GINAN/r48-continuity-bridge-20260912'
spec['baseline_snapshot_commit']='820c2cb09a6e48a155b3fbdc69d100631f88dc6c'
spec['baseline_outputs']['r48fix1']='/mnt/d/GINAN_R20/inputData/outputs/'+previous
spec['checkpoints']['directory']=checkpoint
spec['source_chat']={'project':'OSB','title':'分析解算速度并测试','id':'6aa68697-6fac-83e9-be18-acaf20169431','latest_answer':'9e70982f-4190-4cda-9a00-42d9194d40d0','read_turns':4,'has_more':False,'text_truncated':False,'attachment_contents_retrieved':False,'transcript_sha256':hashlib.sha256((w/'source_conversation.md').read_bytes()).hexdigest()}
spec.pop('repair_parent_commit',None);spec.pop('repair',None)
spec['r49_policy']={'formal_stage':'R49-D','fusion_weight':0.10,'fresh_weight':0.45,'history_weight':0.45,'bridge_search_cap_per_route':4,'fusion_search_cap':4,'partial_shadow_formal':False,'partial_shadow_separate_test':True,'no_statistical_threshold_relaxation':True,'actual_pea_blas_audit':True}
spec['budget']['route_allocation']='0.45/0.45/0.10 fresh/history/fusion, before statistical searches; missing history share goes to fresh'
spec['claims']='Internal integer certification and product continuity only; not independent integer truth or positioning benefit. Compare R47 on 47 common epochs and R48fix1 on 61 common epochs. Separate strategy changes from reference equivalence tests and never sum nested timers.'
(w/'experiment.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2)+'\n')
(w/'config').mkdir(exist_ok=False)
for p in (old/'config').iterdir():
 if p.is_file():
  text=p.read_text().replace(previous,case).replace('/mnt/d/GINAN_R20/r48fix1_checkpoints_20260913',checkpoint).replace('R48-NETWORK','R49-NETWORK')
  (w/'config'/p.name).write_text(text)
for name in ['INPUTS_SHA256SUMS','metrics.py']:
 shutil.copyfile(old/name,w/name)
for oldname,newname in [('smoke.py','smoke.py'),('verify_build.py','verify_build.py'),('launch_r48.py','launch_r49.py'),('run_r48_30m.sh','run_r49_30m.sh')]:
 text=(old/oldname).read_text().replace(str(old),str(w)).replace('r48fix1_work_20260913','r49_work_20260913').replace('r48-continuity-bridge-20260912','r49-condition-domain-20260913').replace('frozen-r48fix1-active-basis-20260913','frozen-r49-condition-domain-20260913').replace('r48fix1_checkpoints_20260913','r49_checkpoints_20260913').replace('r48fix1_checkpoint_smoke_20260913','r49_checkpoint_smoke_20260913').replace(previous,case).replace('frozen_r48.json','frozen_r49.json').replace('run_r48_30m.sh','run_r49_30m.sh').replace('verify_r48.py','verify_r49.py').replace('r48_', 'r49_')
 if oldname=='run_r48_30m.sh':
  text=text.replace('export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1','export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1\nexport ZHANG_R49_FUSION_WEIGHT=0.10\nunset ZHANG_R49_PARTIAL_SHADOW ZHANG_R49_BLAS_INJECT_INVALID')
 if oldname=='smoke.py':
  text=text.replace("'MKL_NUM_THREADS':'1'","'MKL_NUM_THREADS':'1','ZHANG_R49_FUSION_WEIGHT':'0.10'")
 (w/newname).write_text(text)
print(json.dumps({'case':case,'work':str(w),'source':source,'config_count':len(list((w/'config').glob('*.yaml'))),'input_data_copied':False},indent=2))
