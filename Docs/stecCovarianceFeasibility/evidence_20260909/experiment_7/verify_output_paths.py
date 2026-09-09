from pathlib import Path
import csv, json, sys
r=Path('/mnt/d/tec/code-exp7-joint-gate-20260909')
b=Path('/mnt/d/tec/code-exp6-joint-nis-20260909')
mode=sys.argv[1]
decisions=json.loads((r/f'{mode}_gate_verification.json').read_text())
def tow(t):
    h,m,s=map(int,t[11:].split(':')); return 259200+3600*h+60*m+s
passed=set(map(tow,decisions['passed_epochs_GPST']))
rejected=set(map(tow,decisions['rejected_epochs_GPST']))
def states(directory):
    p=next((directory/'ionstec').glob('*.STEC.COV'))
    out={}
    for row in csv.reader(p.open()):
        if row and row[0]=='STATE':
            t=int(float(row[2]))
            if t>334290: break
            out[(t,row[4],row[5],row[6])]=(float(row[8]),float(row[9]))
    return out
actual=states(r/f'replay_{mode}')
baseline=states(b/f'replay_{mode}')
def compare(reference,times):
    selected={k:v for k,v in actual.items() if k[0] in times}
    assert set(selected)=={k for k in reference if k[0] in times}
    errx=max((abs(v[0]-reference[k][0]) for k,v in selected.items()),default=0.)
    errp=max((abs(v[1]-reference[k][1]) for k,v in selected.items()),default=0.)
    assert errx<1e-5 and errp<1e-5, (errx,errp)
    return {'epochs':len(times),'states':len(selected),'max_STEC_error_TECU':errx,'max_variance_error_TECU2':errp}
result={'passed_vs_original_AR':compare(baseline,passed)}
if mode=='10mm':
    floating=states(Path('/mnt/d/tec/code-exp5-sigma-20260909-v2/sigma_10mm_float'))
    result['rejected_vs_independent_FLOAT']=compare(floating,rejected)
(r/f'{mode}_output_path_verification.json').open('x').write(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
