from pathlib import Path
import json,re,collections
W=Path(__file__).resolve().parent
tags={'ZHANG_PRODUCT_INTEGER_LEDGER','ZHANG_PRODUCT_PUBLICATION_AUTHORITY','ZHANG_R47_LEDGER_TRANSACTION','ZHANG_R47_WRITER_ATOMIC_COMMIT','ZHANG_PRODUCT_LATTICE_POST_LEDGER_GRAPH','ZHANG_PRODUCT_LEDGER_CANONICALIZATION','ZHANG_CURRENT_EPOCH_CERTIFIED_LATTICE','ZHANG_R47_SEARCH_ROUTE','R49_FINAL_DOMAIN','ZHANG_PRODUCT_INTEGER_LEDGER_CONFLICT_DETAIL'}
records=[];text=[]
with (W/'final_zero_evidence.txt').open() as f:
 for line in f:
  m=re.match(r'(R49|R48fix1):([^:]+):(\d+): (\w+) time=2024-07-17 (\d\d:\d\d:\d\d)',line)
  if not m or m[4] not in tags:continue
  d={'version':m[1],'trace':m[2],'line':int(m[3]),'tag':m[4],'time':m[5], 'fields':dict(re.findall(r'\b(\w+)=([^\s]+)',line[m.end():]))}
  records.append(d)
  if m[5]>='00:29:00': text.append(line.rstrip())
(W/'final_zero_key_evidence.txt').write_text('\n'.join(text)+'\n')
(W/'final_zero_diagnosis.json').write_text(json.dumps(records,indent=2)+'\n')
for r in records:
 if r['time']<'00:29:00':continue
 fields=r['fields']
 keys=['input_rows','fresh_rows','confirmed_rows','conflicting_rows','status','reason','committed','certified_pair_rank','applied_conditioning_rank','projected_rows','selection_nis','selection_nis_threshold','selection_maximum_perr','current_conditioning_rank']
 if r['tag'] not in {'ZHANG_PRODUCT_INTEGER_LEDGER','ZHANG_PRODUCT_LEDGER_CANONICALIZATION','ZHANG_PRODUCT_LATTICE_POST_LEDGER_GRAPH','ZHANG_R47_WRITER_ATOMIC_COMMIT'}:continue
 print(r['version'],r['time'],r['tag'],json.dumps({k:fields[k] for k in keys if k in fields}))
