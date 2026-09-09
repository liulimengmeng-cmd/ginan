from pathlib import Path
import sys,re,json,numpy as np
sys.path.insert(0,'/mnt/d/tec/ginan-main-code-products/scripts')
import audit_pppar_same_epoch_joint_nis as audit
from audit_pppar_feedback_information import fields
r=Path('/mnt/d/tec/code-exp6-joint-nis-20260909')
mode=sys.argv[1]
original=Path('/mnt/d/tec/code-exp5-sigma-20260909-v2/sigma_10mm_ar') if mode=='10mm' else Path('/mnt/d/tec/code-exp4c-full-20260909')
trace=next(x for x in original.glob('Network-*.trace') if 'smoothed' not in x.name)
snapshots={audit.read_snapshot(p)[0]:p for p in (r/f'snapshots_{mode}').glob('*.snapshot')}
events={t:{'candidate_rows':[],'stec_AR':[]} for t in snapshots}
gpst=None; ar=False
for line in trace.open(errors='replace'):
    if line.startswith('fixAndHoldAmbiguities:'):
        gpst=line.split(': ',1)[1].strip()[:19]
    if '+STATES/AR' in line: ar=True
    if '-STATES/AR' in line: ar=False
    if gpst not in events: continue
    ev=events[gpst]
    if 'PPP_AR DUAL_FREQUENCY_CANDIDATE_ROW ' in line:
        d=fields(line)
        if d['candidate_scope'] in ('SELECTED_SUBSET','FULL_VISIBLE_GROUP'): ev['candidate_rows'].append(d)
    if 'PPP_AR FEEDBACK_SHADOW_SUMMARY ' in line: ev['shadow']=fields(line)
    if ar and line.startswith('*'):
        cols=[x.strip() for x in line.split('\t')]
        if len(cols)>9 and cols[3]=='IONO_STEC' and cols[2][:19]==gpst:
            ev['stec_AR'].append({'receiver':cols[5],'satellite':cols[4],'state':float(cols[7]),'adjust':float(cols[9])})
    if gpst>'2024-07-17 20:51:30': break
reports={}
for t,p in sorted(snapshots.items()):
    _,m,k=audit.read_snapshot(p); result=audit.analyse(m,k); ev=events[t]
    original_nis=float(ev['shadow']['joint_nis'])
    err=abs(original_nis-result['full_nis'])
    assert err<1e-7*max(1,original_nis)
    rows=ev['candidate_rows']; assert len(rows)==m['H'].shape[0]
    expected_v=np.array([-float(d['float_minus_integer']) for d in rows])
    verr=float(np.max(np.abs(expected_v-m['v'].ravel())))
    predicted_sigma=np.sqrt(np.diag(m['H']@m['P']@m['H'].T))
    serr=float(np.max(np.abs(predicted_sigma-np.array([float(d['formal_sigma']) for d in rows]))))
    assert verr<1e-9 and serr<1e-9
    report={'original_joint_nis':original_nis,'snapshot_joint_nis':result['full_nis'],'difference':err,
            'canonical_candidate_to_raw_filter_residual_error_max':verr,'candidate_sigma_vs_snapshot_error_max':serr,
            'rows':len(rows),'subsets':{name:{q:z[q] for q in ('rows','nis','nis_per_row')} for name,z in result['subsets'].items()}}
    out=r/(p.stem+f'_{mode}_joint_verified.json')
    with out.open('x') as f: json.dump(result,f,indent=2)
    if mode=='10mm' and t=='2024-07-17 20:51:00':
        S=m['H']@m['P']@m['H'].T+m['R']; S[9,9]+=m['R'][9,9]*(1e6-1)
        dx=m['P']@m['H'].T@np.linalg.solve(S,m['v'].ravel())
        comparisons=[]
        for row in ev['stec_AR']:
            idx=next(i for i,key in k.items() if key['type']=='IONO_STEC' and key['receiver']==row['receiver'] and key['satellite']==row['satellite'])
            comparisons.append({**row,'prior':float(m['x'][idx,0]),'observed_change':row['state']-float(m['x'][idx,0]),'replayed_change':float(dx[idx]),'replay_error':row['state']-float(m['x'][idx,0])-float(dx[idx])})
        assert len(comparisons)==50
        report['actual_stec_AR_vs_float']=comparisons
        report['stec_replay_error_max']=max(abs(x['replay_error']) for x in comparisons)
        assert report['stec_replay_error_max']<1e-6
    reports[t]=report
(r/f'{mode}_original_replay_verification.json').open('x').write(json.dumps(reports,indent=2))
print(json.dumps({t:{k:v for k,v in x.items() if k!='actual_stec_AR_vs_float'} for t,x in reports.items()},indent=2))
