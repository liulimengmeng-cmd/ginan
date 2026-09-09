from pathlib import Path
import os,sys,subprocess,json,time,resource,hashlib
r=Path('/mnt/d/tec/code-exp6-joint-nis-20260909')
repo=Path('/mnt/d/tec/ginan-main-code-products')
mode=sys.argv[1]
assert mode in ('10mm','0mm')
cfg0=Path('/mnt/d/tec/code-exp5-sigma-20260909-v2/sigma_10mm_ar.yaml') if mode=='10mm' else repo/'Docs/stecCovarianceFeasibility/experiment_4c_2880e_code_rapid_persistent_phase_osb.yaml'
cfg=r/f'replay_{mode}.yaml'
text=cfg0.read_text()
assert text.count('    chunking:')==1
cfg.open('x').write(text.replace('    chunking:','    rts:\n      enable: false\n    chunking:'))
out=r/f'replay_{mode}'; out.mkdir(exist_ok=False)
snaps=r/f'snapshots_{mode}'; snaps.mkdir(exist_ok=True)
times=['2024-07-17 05:13:30','2024-07-17 20:50:00','2024-07-17 20:50:30','2024-07-17 20:51:00','2024-07-17 20:51:30']
env=os.environ.copy(); env.update(LD_LIBRARY_PATH='/home/rx/.local/boost-1.82/lib',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',GINAN_AR_SNAPSHOT_DIRECTORY=str(snaps),GINAN_AR_SNAPSHOT_TIMES='|'+'|'.join(times)+'|')
argv=[str(r/'pea_snapshot'),'-y',str(cfg),'-n','2504','-a','EXP1_DATA_ROOT:/home/rx/GINAN/inputData','-a',f'EXP1_OUTPUT_ROOT:{out}']
def limits():
    os.nice(19)
    resource.setrlimit(resource.RLIMIT_AS,(2200*1024**2,2200*1024**2))
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
start=time.time()
with (r/f'{mode}.stdout.log').open('x') as so,(r/f'{mode}.stderr.log').open('x') as se:
    p=subprocess.Popen(argv,cwd=repo,env=env,stdout=so,stderr=se,preexec_fn=limits)
    (r/f'{mode}.started.json').open('x').write(json.dumps({'pid':p.pid,'argv':argv,'epoch_time_system':'GPST','snapshot_times':times,'changes_from_exp5':'read-only snapshots, RTS output disabled, stop at 2504 epochs','binary_sha256':sha(r/'pea_snapshot'),'config_sha256':sha(cfg),'nice':19,'address_space_limit_MiB':2200,'existing_other_pea_pid':522736,'started_unix':start},indent=2))
    rc=p.wait()
(r/f'{mode}.finished.json').open('x').write(json.dumps({'returncode':rc,'elapsed_seconds':time.time()-start,'binary_sha256':sha(r/'pea_snapshot'),'snapshots':[{'path':str(p),'sha256':sha(p)} for p in snaps.glob('*.snapshot')]},indent=2))
print(mode,rc,flush=True)
raise SystemExit(rc)
