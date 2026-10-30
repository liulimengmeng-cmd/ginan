"""Independent RINEX/SP3/CLK fixture generation; no input modifications.

Transmit-time nine-point SP3 interpolation + Sagnac + relativistic clock.
This is not a replacement for the production ephemeris/antenna model: the
production-adapter test and a controlled PEA replay are separate acceptance gates.
"""
import argparse
import bisect
import hashlib
import json
import math
from pathlib import Path

C = 299792458.0
OMEGA = 7.2921151467e-5


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def read_rinex(path, start, end):
    types, system, position = {}, None, None
    rows = []
    with path.open() as stream:
        for line in stream:
            label = line[60:].strip()
            if label == 'APPROX POSITION XYZ':
                position = [float(x) for x in line[:60].split()]
            if label == 'SYS / # / OBS TYPES':
                if line[0] != ' ':
                    system = line[0]
                    types[system] = []
                types[system].extend(line[7:60].split())
            if label == 'END OF HEADER':
                break
        time = None
        for line in stream:
            if line.startswith('>'):
                f = line.split()
                time = int(f[4]) * 3600 + int(f[5]) * 60 + float(f[6])
                if time > end:
                    break
            elif time is not None and time >= start and line.startswith('G'):
                obs = {}
                for i, code in enumerate(types['G']):
                    field = line[3+16*i:3+16*(i+1)].rstrip().ljust(16)
                    if field[:14].strip():
                        obs[code] = (float(field[:14]), int(field[14]) if field[14].isdigit() else 0)
                if all(obs.get(code, (0, 0))[0] != 0 for code in ['L1C', 'L2W', 'C1C', 'C2W']):
                    rows.append((time, int(line[1:3]), obs))
    return position, rows


def read_sp3(path):
    points = {}
    with path.open() as stream:
        for line in stream:
            if line.startswith('*'):
                f = line.split()
                time = int(f[4])*3600+int(f[5])*60+float(f[6])
            if line.startswith('PG'):
                prn = int(line[2:4])
                points.setdefault(prn, []).append((time, [float(line[a:b])*1000 for a, b in [(4,18),(18,32),(32,46)]]))
    return points


def read_clock(path, end):
    points = {}
    with path.open() as stream:
        for line in stream:
            if not line.startswith('AS G'):
                continue
            f = line.split()
            time = int(f[5])*3600+int(f[6])*60+float(f[7])
            if time > end+300:
                break
            points.setdefault(int(f[1][1:]), []).append((time,float(f[9])))
    return points


def interpolate(points, time):
    near = sorted(sorted(points, key=lambda row: abs(row[0]-time))[:9])
    result = [0.,0.,0.]
    for i, (ti, xyz) in enumerate(near):
        w = math.prod((time-tj)/(ti-tj) for j, (tj, _) in enumerate(near) if i != j)
        for k in range(3):
            result[k] += w*xyz[k]
    return result


def clock_value(points, time):
    i = min(max(1,bisect.bisect_left([p[0] for p in points], time)), len(points)-1)
    t0,y0 = points[i-1];t1,y1 = points[i]
    return y0+(y1-y0)*(time-t0)/(t1-t0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    rinex = args.input_root/'data/MAL200KEN_R_20241990000_01D_30S_MO.rnx'
    sp3 = args.input_root/'products/WUM0MGXFIN_20241990000_01D_05M_ORB.SP3'
    clk = args.input_root/'products/WUM0MGXFIN_20241990000_01D_30S_CLK.CLK'
    position, rows = read_rinex(rinex, 6600, 7800)
    ephemeris = read_sp3(sp3); clocks = read_clock(clk,7800)
    horizontal = math.hypot(*position[:2]); lat=math.atan2(position[2],horizontal*(1-6.69437999014e-3))
    for _ in range(8):
        radius=6378137/math.sqrt(1-6.69437999014e-3*math.sin(lat)**2)
        lat=math.atan2(position[2]+6.69437999014e-3*radius*math.sin(lat),horizontal)
    lon=math.atan2(position[1],position[0])
    up=[math.cos(lat)*math.cos(lon),math.cos(lat)*math.sin(lon),math.sin(lat)]
    count=0
    with (args.output/'mal2_raw_fixture.txt').open('w') as stream:
        stream.write('# time sat range_minus_satclock elevation L1C L2W L5Q C1C C2W phaseVar codeVar LLI1 LLI2 LLI5 D1C D2W\n')
        for time,sat,obs in rows:
            if sat not in ephemeris or sat not in clocks:
                continue
            transmit=time-obs['C1C'][0]/C
            clock=clock_value(clocks[sat],transmit)
            transmit-=clock
            xyz=interpolate(ephemeris[sat],transmit)
            lo=interpolate(ephemeris[sat],transmit-.01);hi=interpolate(ephemeris[sat],transmit+.01)
            velocity=[(b-a)/.02 for a,b in zip(lo,hi)]
            # GINAN adjustRelativity convention: subtract 2*r.v/c^2.
            clock-=2*sum(a*b for a,b in zip(xyz,velocity))/C**2
            delta=[a-b for a,b in zip(xyz,position)];distance=math.sqrt(sum(a*a for a in delta))
            elevation=math.degrees(math.asin(sum(a*b for a,b in zip(delta,up))/distance))
            if elevation<10:
                continue
            distance+=OMEGA*(xyz[0]*position[1]-xyz[1]*position[0])/C
            factor=math.sin(math.radians(elevation))
            values=[time,sat,distance-C*clock,elevation,obs['L1C'][0],obs['L2W'][0],obs.get('L5Q',(0,0))[0],
                    obs['C1C'][0],obs['C2W'][0],(.002/factor)**2,(.4/factor)**2,
                    obs['L1C'][1],obs['L2W'][1],obs.get('L5Q',(0,0))[1],
                    obs.get('D1C',(0,0))[0],obs.get('D2W',(0,0))[0]]
            stream.write(' '.join(str(value) for value in values)+'\n');count+=1
    manifest={'input_hashes':{str(p):sha(p) for p in [rinex,sp3,clk]},'rows':count,'position':position,
              'scope':'MAL2 GPS all observable links 01:50-02:10, independent geometry fixture; not full PEA replay',
              'production_replay_required':True}
    (args.output/'fixture_manifest.json').write_text(json.dumps(manifest,indent=2))
    print(json.dumps(manifest,indent=2))


if __name__=='__main__':main()
