from pathlib import Path
import numpy as np, scipy.linalg as la,json
w=Path('/mnt/c/Users/rx/Documents/GINAN/r48_user_model_work_20260914');results=[]
def read(p):
 d={}
 with p.open('rb') as f:
  assert f.readline()==b'USER_MEASUREMENT_V1\n';label=f.readline().decode().strip()
  while line:=f.readline():
   name,n,m=line.decode().split();n=int(n);m=int(m);d[name]=np.frombuffer(f.read(n*m*8),dtype='<f8').reshape((n,m),order='F').copy();assert f.read(1)==b'\n'
 return d,label
for station,model in [('MARS','UDUC'),('BREW','IF')]:
 for joseph in [0,1]:
  case=f'r48u_20260914_dump_{station}_{model}_{joseph}';root=Path('/mnt/d/GINAN_R20/inputData/outputs')/case/'measurement_dumps'
  files=sorted(root.glob('*.bin'),key=lambda p:int(p.stem.split('_')[-1]))
  if not files:continue
  d,label=read(files[-1]);p=d['P'];h=d['H'];r=d['R'];v=d['v'];n=len(p)
  keep=np.flatnonzero(np.diag(p)>0);sub=p[np.ix_(keep,keep)];scale=np.sqrt(np.diag(sub));corr=sub/np.outer(scale,scale)
  L=np.zeros((n,len(keep)));L[keep,:]=scale[:,None]*la.cholesky((corr+corr.T)/2,lower=True)
  lr=la.cholesky((r+r.T)/2,lower=True);A=la.solve_triangular(lr,h@L,lower=True)
  aug=np.vstack([np.eye(len(keep)),A]);U,T=la.qr(aug,mode='economic');b=np.vstack([np.zeros((len(keep),1)),la.solve_triangular(lr,v,lower=True)])
  z=la.solve_triangular(T,U.T@b);rootcov=la.solve_triangular(T.T,L.T,lower=True).T;qrP=rootcov@rootcov.T;qrx=d['x']+L@z
  S=h@p@h.T+r;ss=np.sqrt(np.diag(S));Sc=S/np.outer(ss,ss);K=la.solve(Sc,(h@p)/ss[:,None],assume_a='sym').T/ss[None,:]
  I=np.eye(n)-K@h;jp=I@p@I.T+K@r@K.T
  denom=max(1,la.norm(jp));res={'case':case,'epoch':label,'n':n,'m':len(r),'condition_S':float(np.linalg.cond(S)),'KF_vs_Joseph':float(la.norm(d['Pp']-jp)/denom),'KF_vs_QR':float(la.norm(d['Pp']-qrP)/max(1,la.norm(qrP))),'Joseph_vs_QR':float(la.norm(jp-qrP)/denom),'KF_mean_vs_QR':float(la.norm(d['xp']-qrx)/max(1,la.norm(qrx))),'min_QR_eigenvalue':float(la.eigvalsh(qrP)[0])}
  results.append(res);print(json.dumps(res),flush=True)
(w/'independent_posterior_audit.json').write_text(json.dumps(results,indent=2)+'\n')
