"""Run production detector core on raw-derived fixtures and report, not hide, misses."""
import argparse
import csv
import json
import itertools
import os
import subprocess
from pathlib import Path


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--binary',type=Path,required=True)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--grid',action='store_true')
    parser.add_argument('--dual',action='store_true')
    parser.add_argument('--without-doppler',action='store_true')
    args=parser.parse_args()
    cases=[(0,0,0),(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(1,1,0),(-1,-1,0),
           (1,1,1),(-1,-1,-1),(9,7,0),(-9,-7,0),(77,60,0),(60,77,0),(0,0,1),(0,0,-1)]
    fixture=args.root/'mal2_raw_fixture.txt'
    suffix=''
    if args.grid:
        cases=list(dict.fromkeys(cases+list(itertools.product(range(-2,3),repeat=3))))
        suffix='_grid'
    if args.without_doppler:suffix+='_no_doppler'
    if args.dual:
        cases=[case for case in cases if case[2]==0]
        suffix+='_dual'
        fixture=args.root/'mal2_raw_fixture_dual.txt'
        with fixture.open('w') as stream:
            for line in (args.root/'mal2_raw_fixture.txt').read_text().splitlines():
                if not line or line.startswith('#'):
                    stream.write(line+'\n');continue
                fields=line.split();fields[6]='0';fields[13]='0';stream.write(' '.join(fields)+'\n')
    summary=[]
    for n1,n2,n3 in cases:
        name=f'inject{suffix}_{n1}_{n2}_{n3}'
        output=args.root/(name+'.csv')
        with output.open('w') as stream:
            subprocess.run([str(args.binary),'--replay',str(fixture),
                            str(n1),str(n2),str(n3),'7260'],stdout=stream,check=True,
                            env={**os.environ,**({'ZHANG_SLIP_TEST_NO_DOPPLER':'1'} if args.without_doppler else {})})
        with output.open() as stream:
            rows=[r for r in csv.DictReader(stream) if r['satellite']=='25']
        affected=[r for r in rows if float(r['time'])>=7260]
        detected=next((r for r in affected if r['use_phase']=='0' or r['break_arc']=='1'),None)
        confirmed=next((r for r in affected if r['action']=='INTEGER_STEP'),None)
        unknown=next((r for r in affected if r['action']=='CONTINUITY_UNKNOWN_PHASE_STEP'),None)
        record={'injection_cycles':[n1,n2,n3],'detected':detected is not None,
                'first_detection_seconds':float(detected['time'])-7260 if detected else None,
                'confirmed':confirmed is not None,'initial_action':affected[0]['action'],
                'continuity_unknown_retirement':unknown is not None,
                'initial_phase_score':float(affected[0]['phase_score']),
                'confirmed_integer':[int(confirmed[k]) for k in ['n1','n2','n3']] if confirmed else None,
                'quarantines':sum(r['use_phase']=='0' for r in affected)}
        summary.append(record)
    result={'cases':summary,'all_injected_cases_detected':all(r['detected'] for r in summary[1:]),
            'doppler_intentionally_unavailable':args.without_doppler,
            'MAL2_G25_no_injection_quarantines':summary[0]['quarantines'],
            'restart_authorized_by_tests':False,
            'scope':'detector core + independent geometry only; production adapter/checkpoint tests still required'}
    result['wrong_integer_certificates']=[r for r in summary[1:]
        if r['confirmed'] and r['confirmed_integer']!=r['injection_cycles']]
    (args.root/('injection_summary'+suffix+'.json')).write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
    return 0 if (result['all_injected_cases_detected'] and
                 result['MAL2_G25_no_injection_quarantines']==0 and
                 not result['wrong_integer_certificates']) else 1


if __name__=='__main__':raise SystemExit(main())
