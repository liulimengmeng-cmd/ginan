"""Strict committed-product/decision comparison; no re-alignment or fitting.
OLD_TRACE NEW_TRACE OLD_BUNDLE NEW_BUNDLE OUTPUT.json
"""
from pathlib import Path
import csv
import json
import math
import re
import sys

old_trace,new_trace,old_bundle,new_bundle,output=map(Path,sys.argv[1:6])
def decisions(path):
    result=[]
    for line in path.open():
        if line.startswith('ZHANG_R51_BLOCK_ACCEPT '):
            fields=dict(re.findall(r'(\w+)=(\S+)',line))
            result.append({k:fields[k] for k in ('stage','scan','accepted_rank','total','ils_calls','ratio_threshold')})
    return result
def table(directory,name,keys):
    with (directory/name).open() as stream:rows=list(csv.DictReader(stream))
    result={tuple(r[k] for k in keys):r for r in rows}
    assert len(result)==len(rows),'duplicate keys'
    return result

keys=('gpst_seconds','solution','satellite','observable')
a,b=(table(p,'products.csv',keys) for p in (old_bundle,new_bundle))
assert a.keys()==b.keys(),'product membership changed'
numeric=('clock_m','phase_m','correction_m','clock_sigma_m','phase_sigma_m',
         'correction_sigma_m','clock_phase_covariance_m2')
maxima={}
for field in numeric:
    differences=[]
    for key in a:
        x,y=float(a[key][field]),float(b[key][field])
        assert math.isfinite(x)==math.isfinite(y),f'finite mask changed: {key} {field}'
        if math.isfinite(x):differences.append(abs(x-y))
    maxima[field]=max(differences,default=0.)
flags=('ar_valid','pppar_usable','ppp_usable','product_state','discontinuity_counter',
       'discontinuity','public_phase_segment_id','integer_alignment_generation',
       'real_gauge_generation','same_segment_continuous','requires_user_ambiguity_reset',
       'integer_structure_valid','integer_datum_continuous','integer_precision_valid',
       'integer_valid','cycle_closure_valid','dual_frequency_ar_valid','joint_covariance_complete',
       'continuity_valid','numeric_valid','branch_valid','float_chart_transport_valid',
       'integer_transport_valid','integer_transport_shift_cycles','integer_component_rank',
       'certified_relation_count','redundant_relation_count','integer_shift_cycles','fractional_shift_cycles')
flag_changes={f:sum(a[k][f]!=b[k][f] for k in a) for f in flags}
all_field_changes={f:sum(a[k][f]!=b[k][f] for k in a) for f in next(iter(a.values()))}
covkeys=('gpst_seconds','solution','row_satellite','row_parameter','row_observable',
         'column_satellite','column_parameter','column_observable')
ac,bc=(table(p,'covariance.csv',covkeys) for p in (old_bundle,new_bundle))
assert ac.keys()==bc.keys(),'covariance membership changed'
cov_error=0.
for key in ac:
    x,y=float(ac[key]['covariance_m2']),float(bc[key]['covariance_m2'])
    assert math.isfinite(x)==math.isfinite(y),'covariance finite mask changed'
    if math.isfinite(x):cov_error=max(cov_error,abs(x-y))
da,db=decisions(old_trace),decisions(new_trace)
result=dict(old_bundle=str(old_bundle),new_bundle=str(new_bundle),
    old_trace=str(old_trace),new_trace=str(new_trace),product_rows=len(a),covariance_entries=len(ac),
    maximum_product_differences=maxima,maximum_covariance_difference=cov_error,
    publication_and_integer_flag_changes=flag_changes,all_field_changes=all_field_changes,
    accepted_block_decisions_equal=da==db,accepted_block_count=len(db),
    passed=da==db and not any(flag_changes.values()) and
        max(maxima.values())<1e-7 and cov_error<1e-9)
output.write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
assert result['passed'],'real epoch comparison failed; retain evidence'
