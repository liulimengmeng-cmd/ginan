"""Stream every old/new mean and covariance entry; do not fit offsets.
OLD NEW OUTPUT [suite]. Default is a single-dimension production replay.
"""
from pathlib import Path
import json
import struct
import sys
import numpy as np

old,new,output=map(Path,sys.argv[1:4]);suite=len(sys.argv)>4
def integer(stream):
    data=stream.read(4)
    if not data:return None
    assert len(data)==4
    return struct.unpack('i',data)[0]
def record(stream,n):
    valid=integer(stream)
    if valid is None:return None
    rank,det,length=(integer(stream) for _ in range(3))
    reason=stream.read(length).decode()
    values=np.frombuffer(stream.read((n+n*n)*8),dtype=np.float64) if valid else np.empty(0)
    assert not valid or (len(values)==n+n*n and np.isfinite(values).all())
    return (valid,rank,det,reason),values

max_mean=max_cov=0.;compared=entries=0;status_changes=[]
with old.open('rb') as a,new.open('rb') as b:
    if not suite:
        n=integer(a);assert n==integer(b)
    while True:
        if suite:
            n=integer(a);assert n==integer(b)
            if n is None:break
        left,right=record(a,n),record(b,n)
        assert (left is None)==(right is None),'different record counts'
        if left is None:break
        compared+=1
        if left[0]!=right[0]:status_changes.append(dict(case=compared,old=left[0],new=right[0]))
        if left[0][0] and right[0][0]:
            error=np.abs(left[1]-right[1]);entries+=len(error)
            max_mean=max(max_mean,float(error[:n].max(initial=0)))
            max_cov=max(max_cov,float(error[n:].max(initial=0)))
result=dict(old=str(old),new=str(new),records=compared,compared_entries=entries,
    status_changes=status_changes,maximum_mean_difference=max_mean,
    maximum_covariance_difference=max_cov,tolerance=2e-10,
    passed=not status_changes and max(max_mean,max_cov)<2e-10)
output.write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
assert result['passed'],'conditioner comparison failed; preserve failure evidence'
