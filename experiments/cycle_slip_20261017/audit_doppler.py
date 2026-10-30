"""Independent current-epoch phase/Doppler audit; all inputs remain read-only."""
import argparse
import json
import math
from pathlib import Path
import statistics
from prepare_replay import read_rinex, read_sp3, interpolate, C, OMEGA, sha


def geometry(points,position,time):
    travel=.075
    for _ in range(4):
        xyz=interpolate(points,time-travel)
        distance=math.sqrt(sum((a-b)**2 for a,b in zip(xyz,position)))
        distance+=OMEGA*(xyz[0]*position[1]-xyz[1]*position[0])/C
        travel=distance/C
    return distance


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    rinex=args.input_root/'data/MAL200KEN_R_20241990000_01D_30S_MO.rnx'
    sp3=args.input_root/'products/WUM0MGXFIN_20241990000_01D_05M_ORB.SP3'
    position,rows=read_rinex(rinex,6000,8100)
    ephemeris=read_sp3(sp3)
    frequencies=[('L1C','D1C',C/1575420000),('L2W','D2W',C/1227600000),('L5Q','D5Q',C/1176450000)]
    previous={};epochs={};missing=0
    for time,sat,obs in rows:
        if sat not in ephemeris:continue
        distance=geometry(ephemeris[sat],position,time)
        rate=(geometry(ephemeris[sat],position,time+.1)-geometry(ephemeris[sat],position,time-.1))/.2
        prior=previous.get(sat);previous[sat]=(time,obs,distance,rate)
        if not prior or time-prior[0]!=30:continue
        pt,po,pr,pv=prior;dt=time-pt
        trap=distance-pr-.5*dt*(rate+pv)
        residual=[]
        for phase,doppler,wavelength in frequencies:
            if not all(key in data and data[key][0]!=0 for data in [obs,po] for key in [phase,doppler]):
                residual.append(None);missing+=1;continue
            if obs[phase][1]&3 or po[phase][1]&3:residual.append(None);continue
            residual.append(wavelength*(obs[phase][0]-po[phase][0]+.5*dt*(obs[doppler][0]+po[doppler][0]))-trap)
        epochs.setdefault(time,[]).append({'satellite':sat,'trap_correction_m':trap,'raw_error_m':residual})
    all_errors=[[] for _ in frequencies];target=[];gf_errors=[]
    for time,entries in epochs.items():
        for entry in entries:
            e=entry['raw_error_m']
            if e[0] is not None and e[1] is not None:gf_errors.append(e[0]-e[1])
            errors=[]
            for f,value in enumerate(entry['raw_error_m']):
                others=[r['raw_error_m'][f] for r in entries if r is not entry and r['raw_error_m'][f] is not None]
                error=value-statistics.median(others) if value is not None and len(others)>=3 else None
                errors.append(error)
                if error is not None:all_errors[f].append(error)
            if entry['satellite']==25:target.append({'time':time,**entry,'clock_removed_error_m':errors})
    def describe(values):
        values=sorted(values);median=statistics.median(values)
        return {'n':len(values),'median_m':median,'mad_sigma_m':1.4826*statistics.median(abs(v-median) for v in values),
            'max_abs_m':max(abs(v) for v in values),'p95_abs_m':sorted(abs(v) for v in values)[int(.95*(len(values)-1))]}
    result={'scope':'independent raw RINEX/SP3 diagnostic, not production KF/graph replay',
        'rinex_sha256':sha(rinex),'sp3_sha256':sha(sp3),'missing_doppler_pairs':missing,
        'clock_removed_summary':[describe(v) for v in all_errors],
        'geometry_free_doppler_summary':describe(gf_errors),'MAL2_G25':target}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2))
    print(json.dumps({'clock_removed_summary':result['clock_removed_summary'],
        'geometry_free_doppler_summary':result['geometry_free_doppler_summary'],
        'MAL2_G25_event':[v for v in target if 7200<=v['time']<=7410]},indent=2))


if __name__=='__main__':main()
