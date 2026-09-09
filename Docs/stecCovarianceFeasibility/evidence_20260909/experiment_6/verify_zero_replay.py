from pathlib import Path
import sys,re,json,numpy as np
sys.path.insert(0,'/mnt/d/tec/ginan-main-code-products/scripts')
import audit_pppar_same_epoch_joint_nis as audit
from audit_pppar_feedback_information import fields
r=Path('/mnt/d/tec/code-exp6-joint-nis-20260909')
snapshots={audit.read_snapshot(p)[0]:p for p in (r/'snapshots_0mm').glob('*.snapshot')}
def read_events(root):
    p=next(x for x in root.glob('Network-*.trace') if 'smoothed' not in x.name)
    result={t:{'rows':[]} for t in snapshots}; t=None
    for line in p.open(errors='replace'):
        if line.startswith('fixAndHoldAmbiguities:'): t=line.split(': ',1)[1].strip()[:19]
        if t not in result: continue
        if 'PPP_AR DUAL_FREQUENCY_CANDIDATE_ROW ' in line:
            d=fields(line)
            if d['candidate_scope'] in ('SELECTED_SUBSET','FULL_VISIBLE_GROUP'): result[t]['rows'].append(d)
        if 'PPP_AR FEEDBACK_SHADOW_SUMMARY ' in line: result[t]['nis']=float(fields(line)['joint_nis'])
        if t>'2024-07-17 20:51:30': break
    return result
historical=read_events(Path('/mnt/d/tec/code-exp4c-full-20260909'))
replay=read_events(r/'replay_0mm')
reports={}
for t,p in snapshots.items():
    _,m,k=audit.read_snapshot(p); a=audit.analyse(m,k)
    old=historical[t]; new=replay[t]
    err=abs(a['full_nis']-new['nis']); assert err<1e-7
    assert len(old['rows'])==len(new['rows'])==len(a['rows'])
    same_integers=all((o['receiver'],o['family'],o['rhs'],o['reference'])==(n['receiver'],n['family'],n['rhs'],n['reference']) for o,n in zip(old['rows'],new['rows']))
    assert same_integers
    verr=float(np.max(np.abs(m['v'].ravel()+np.array([float(d['float_minus_integer']) for d in new['rows']]))))
    serr=float(np.max(np.abs(np.sqrt(np.diag(m['H']@m['P']@m['H'].T))-np.array([float(d['formal_sigma']) for d in new['rows']]))))
    assert verr<1e-9 and serr<1e-9
    (r/(p.stem+'_0mm_joint_verified.json')).open('x').write(json.dumps(a,indent=2))
    reports[t]={'snapshot_vs_own_replay_nis_error':err,'candidate_residual_error':verr,'candidate_sigma_error':serr,
        'same_historical_integer_candidates':same_integers,'historical_nis':old['nis'],'replay_nis':new['nis'],
        'historical_relative_nis_difference':abs(new['nis']-old['nis'])/old['nis'],
        'historical_float_residual_difference_max':max(abs(float(o['float_minus_integer'])-float(n['float_minus_integer'])) for o,n in zip(old['rows'],new['rows'])),
        'rows':len(a['rows']),'full_nis_per_row':a['full_nis_per_row'],
        'cross_station_covariance_removed_nis_per_row':a['cross_station_covariance_removed_nis']/len(a['rows'])}
    # Do not claim bitwise historical reproduction; quantify the small drift.
    assert reports[t]['historical_relative_nis_difference']<1e-5
(r/'0mm_original_replay_verification.json').open('x').write(json.dumps(reports,indent=2))
print(json.dumps(reports,indent=2))
