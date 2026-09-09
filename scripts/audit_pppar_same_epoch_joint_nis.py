#!/usr/bin/env python3
"""Same-prior marginal/conditional constraint tests from opt-in PEA snapshots.

All row removals and covariance ablations are OFFLINE diagnostics; they neither
select an operational fixed solution nor certify an integer. No covariance floor
is added. Exact posterior filtering with callbacks is not replayed here.
"""
import argparse
from itertools import combinations
import json
from pathlib import Path
import hashlib
import numpy as np


def read_snapshot(path):
    with Path(path).open() as f:
        assert f.readline().strip()=='GINAN_AR_SNAPSHOT_V1'
        epoch=f.readline().strip().removeprefix('GPST ')
        matrices={}
        for name in ('x','P','H','v','R','rhs'):
            label,n,m=f.readline().split(); n,m=int(n),int(m)
            assert label==name
            matrices[name]=np.array([[float(v) for v in f.readline().split()] for _ in range(n)])
            assert matrices[name].shape==(n,m)
        label,n=f.readline().split(); assert label=='KEYS'
        keys={}
        for _ in range(int(n)):
            idx,typ,rec,sat,num=f.readline().split()
            keys[int(idx)]={'type':typ,'receiver':rec,'satellite':sat,'signal_number':int(num)}
        assert not f.read().strip()
    return epoch,matrices,keys


def nis(v,S):
    if len(v)==0: return 0.0
    np.linalg.cholesky(S)
    return float(v@np.linalg.solve(S,v))


def analyse(matrices,keys):
    x,P,H,v,R,rhs=(matrices[k] for k in ('x','P','H','v','R','rhs'))
    x=x.ravel(); v=v.ravel(); rhs=rhs.ravel()
    identity_error=float(np.max(np.abs(v-(rhs-H@x))))
    assert identity_error<1e-9
    S0=H@P@H.T+R
    symmetry_error=float(np.max(np.abs(S0-S0.T)))
    assert symmetry_error<1e-9
    S=(S0+S0.T)/2
    labels=[]; rows=[]
    for i,row in enumerate(H):
        support=np.flatnonzero(row)
        receivers={keys[int(j)]['receiver'] for j in support}
        assert len(receivers)==1 and '-' not in receivers
        labels.append(next(iter(receivers)))
        rows.append({'row_zero_based':i,'receiver':labels[-1],'rhs':float(rhs[i]),
            'float_value':float(row@x),'innovation':float(v[i]),'marginal_sigma':float(np.sqrt(S[i,i])),
            'terms':[{'coefficient':float(row[j]),'state_index':int(j),**keys[int(j)]} for j in support]})
    receivers=sorted(set(labels)); labels=np.array(labels)
    def selected(idx,noise=None):
        idx=np.array(idx,dtype=int)
        si=S[np.ix_(idx,idx)].copy()
        if noise is not None: si+=noise
        vi=v[idx]
        score=nis(vi,si)
        delta=P@H[idx].T@np.linalg.solve(si,vi)
        coord={}; stec={}
        for rec in sorted({k['receiver'] for k in keys.values()}-{'-'}):
            ci=[i for i,k in keys.items() if k['type']=='REC_POS' and k['receiver']==rec]
            ii=[i for i,k in keys.items() if k['type']=='IONO_STEC' and k['receiver']==rec]
            if ci: coord[rec]=float(np.linalg.norm(delta[ci]))
            if ii: stec[rec]={'norm_TECU':float(np.linalg.norm(delta[ii])),
                             'max_abs_TECU':float(np.max(np.abs(delta[ii])))}
        return {'rows':len(idx),'row_indices':idx.tolist(),'nis':score,'nis_per_row':score/len(idx),
                'coordinate_update_norm_m_by_receiver':coord,'stec_update_by_receiver':stec}
    subsets={}
    for n in range(1,len(receivers)+1):
        for rs in combinations(receivers,n):
            idx=np.flatnonzero(np.isin(labels,rs))
            subsets['+'.join(rs)]=selected(idx)
    all_idx=np.arange(len(v)); full=nis(v,S)
    conditional={}
    for rec in receivers:
        g=np.flatnonzero(labels==rec); o=np.flatnonzero(labels!=rec)
        if not len(o): continue
        So=S[np.ix_(o,o)]; Sgo=S[np.ix_(g,o)]
        vc=v[g]-Sgo@np.linalg.solve(So,v[o])
        Sc=S[np.ix_(g,g)]-Sgo@np.linalg.solve(So,Sgo.T)
        score=nis(vc,(Sc+Sc.T)/2)
        marginal_other=nis(v[o],So)
        assert abs(score+marginal_other-full)<1e-6*max(1,full)
        conditional[rec]={'conditional_nis':score,'conditional_nis_per_row':score/len(g),
                          'other_stations_nis':marginal_other,'conditional_innovation':vc.tolist(),
                          'conditional_sigma':np.sqrt(np.diag(Sc)).tolist()}
    eig,U=np.linalg.eigh(S); projected=U.T@v; terms=projected**2/eig
    modes=[]
    for j in np.argsort(terms)[::-1][:5]:
        modes.append({'eigenvalue':float(eig[j]),'projected_innovation':float(projected[j]),
            'nis_contribution':float(terms[j]),'row_loadings':U[:,j].tolist(),
            'loading_energy_by_receiver':{rec:float(np.sum(U[labels==rec,j]**2)) for rec in receivers}})
    # Dropping cross-station covariance is explanatory only, never an estimator fix.
    blockS=S.copy(); blockS[labels[:,None]!=labels[None,:]]=0
    deletion=[]
    for i in range(len(v)):
        z=selected(np.delete(all_idx,i)); z.update(removed_row=i,removed_receiver=str(labels[i]))
        deletion.append(z)
    deweight=[]
    for i in range(len(v)):
        extra=np.zeros_like(R); extra[i,i]=R[i,i]*(1e6-1)
        z=selected(all_idx,extra); z.update(deweighted_row=i,receiver=str(labels[i]))
        deweight.append(z)
    return {'identity_error_max':identity_error,'S_symmetry_error_max':symmetry_error,
        'S_min_eigenvalue':float(eig.min()),'S_condition_number':float(eig.max()/eig.min()),
        'full_nis':full,'full_nis_per_row':full/len(v),'rows':rows,'subsets':subsets,
        'conditional_on_other_stations':conditional,
        'cross_station_covariance_removed_nis':nis(v,blockS),
        'dominant_covariance_modes':modes,
        'leave_one_row_out':sorted(deletion,key=lambda z:z['nis']),
        'one_row_variance_times_1e6':sorted(deweight,key=lambda z:z['nis']),
        'claim_limits':['Diagnostic covariance ablation is not a proposed covariance model.',
                        'Predicted state updates are linear, before posterior callbacks.',
                        'Neither NIS nor lower coordinate shift certifies integer truth or STEC accuracy.']}


def main():
    p=argparse.ArgumentParser(); p.add_argument('snapshot',type=Path); p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); epoch,m,k=read_snapshot(a.snapshot)
    result=analyse(m,k); result.update(gpst=epoch,snapshot=str(a.snapshot),
        snapshot_sha256=hashlib.sha256(a.snapshot.read_bytes()).hexdigest())
    with a.output.open('x') as f: json.dump(result,f,indent=2,allow_nan=False)
    print(json.dumps({'gpst':epoch,'full_nis_per_row':result['full_nis_per_row'],
                      'subsets':{k:(v['rows'],v['nis'],v['nis_per_row']) for k,v in result['subsets'].items()}},indent=2))

if __name__=='__main__': main()
