from pathlib import Path
import subprocess, shlex, json, hashlib, resource, os, time
r=Path('/mnt/d/tec/code-exp6-joint-nis-20260909')
repo=Path('/mnt/d/tec/ginan-main-code-products')
b=Path('/mnt/d/tec/ginan-main-full-rank-ar/build-wsl/cpp')
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
flags={line.split(' = ',1)[0]:line.split(' = ',1)[1] for line in (b/'CMakeFiles/pea.dir/flags.make').read_text().splitlines() if ' = ' in line}
cmd=['/usr/bin/c++']+shlex.split(flags['CXX_DEFINES']+' '+flags['CXX_INCLUDES']+' '+flags['CXX_FLAGS'])
cmd+=['--param','ggc-min-expand=5','--param','ggc-min-heapsize=16384','-o',str(r/'build/ppp_ambres.cpp.o'),'-c',str(repo/'src/cpp/pea/ppp_ambres.cpp')]
def limits():
    os.nice(19)
    resource.setrlimit(resource.RLIMIT_AS,(1400*1024**2,1400*1024**2))
record={'compile_argv':cmd,'baseline_binary_sha256':sha('/mnt/d/tec/ginan-main-full-rank-ar/bin/pea'), 'instrumented_source_sha256':sha(repo/'src/cpp/pea/ppp_ambres.cpp'),'address_space_cap_MiB':1400,'nice':19}
t=time.monotonic()
with (r/'build/compile.log').open('x') as log:
    rc=subprocess.run(cmd,stdout=log,stderr=log,preexec_fn=limits).returncode
record['compile_returncode']=rc
if rc==0:
    link=shlex.split((b/'CMakeFiles/pea.dir/link.txt').read_text())
    record['baseline_link_inputs']=[{'path':str((b/x).resolve()),'sha256':sha(b/x)} for x in link if x.endswith('.o')]
    link[link.index('CMakeFiles/pea.dir/pea/ppp_ambres.cpp.o')]=str(r/'build/ppp_ambres.cpp.o')
    link[link.index('-o')+1]=str(r/'pea_snapshot')
    record['link_argv']=link
    with (r/'build/link.log').open('x') as log:
        rc=subprocess.run(link,cwd=b,stdout=log,stderr=log,preexec_fn=limits).returncode
    record['link_returncode']=rc
    if rc==0: record['diagnostic_binary_sha256']=sha(r/'pea_snapshot')
record['elapsed_seconds']=time.monotonic()-t
(r/'build_receipt.json').write_text(json.dumps(record,indent=2))
print('build return code',rc,flush=True)
raise SystemExit(rc)
