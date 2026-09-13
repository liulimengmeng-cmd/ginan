from pathlib import Path
import json,re,collections
from metrics import quality,metrics
W=Path(__file__).resolve().parent
S=json.loads((W/'experiment.json').read_text())
roots={'R47':Path(S['baseline_outputs']['r47']),'R48fix1':Path(S['baseline_outputs']['r48fix1']),'R49':Path('/mnt/d/GINAN_R20/inputData/outputs')/S['case']}
rows={name:quality(root/'zhang_internal_products.csv') for name,root in roots.items()}
indices={name:{r['epoch']:r for r in seq} for name,seq in rows.items()}
lines=['# R49 final-zero epoch comparison','','Strict dual-frequency AR: PRODUCT_FIXED, AR_VALID, both L1C/L2W valid and same nonempty component. Missing R47 epochs are NA, not zero.','','|GPST|R47 AR / largest / components|R48fix1 AR / largest / components|R49 AR / largest / components|','|---|---:|---:|---:|']
for epoch in sorted(set().union(*(set(v) for v in indices.values()))):
    cells=[]
    for name,index in indices.items():
        r=index.get(epoch)
        cells.append(f"{r['strict_ar']} / {r['largest_component']} / {len(r['components'])}" if r else 'NA')
    r=next(v[epoch] for v in indices.values() if epoch in v)
    lines.append('|'+r['time']+'|'+'|'.join(cells)+'|')
(W/'R49_final_zero_epoch_comparison.md').write_text('\n'.join(lines)+'\n')
evidence=[]
for name,root in roots.items():
    for p in root.glob('Network*.TRACE'):
        with p.open(errors='replace') as f:
            for lineno,line in enumerate(f,1):
                if not line.startswith(('ZHANG_','R48_','R49_')):continue
                m=re.search(r'time=2024-07-17 (\d\d:\d\d:\d\d)',line)
                if m and m[1]>='00:25:00':evidence.append(f'{name}:{p.name}:{lineno}: '+line.rstrip())
(W/'final_zero_evidence.txt').write_text('\n'.join(evidence)+'\n')
print(json.dumps({'rows':{n:len(r) for n,r in rows.items()},'evidence_lines':len(evidence),'metrics':{n:metrics(r) for n,r in rows.items()}},indent=2))
