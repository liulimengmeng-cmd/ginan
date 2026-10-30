"""Inventory the configured network's actual GPS Doppler fields, read-only."""
import argparse
import json
from pathlib import Path
import re


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--input-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    records=[]
    for name in re.findall(r'^\s+- (\S+\.rnx)$',args.config.read_text(),re.M):
        path=args.input_root/name;types=[];system=None
        with path.open() as stream:
            for line in stream:
                label=line[60:].strip()
                if label=='SYS / # / OBS TYPES':
                    if line[0]!=' ':system=line[0]
                    if system=='G':types+=line[7:60].split()
                if label=='END OF HEADER':break
        records.append({'station':name[:4],'path':str(path),'gps_types':types,
            'selected_pair_doppler_declared':all(k in types for k in ['D1C','D2W'])})
    result={'scope':'header declarations only, not per-epoch nonzero availability',
        'stations':len(records),'selected_pair_doppler_declared':sum(r['selected_pair_doppler_declared'] for r in records),
        'missing_stations':[r['station'] for r in records if not r['selected_pair_doppler_declared']],
        'records':records}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k!='records'}))


if __name__=='__main__':main()
