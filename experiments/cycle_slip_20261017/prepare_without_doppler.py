"""Make a read-only-source MAL2 QC fixture with GPS Doppler fields unavailable.

Phase/code/LLI fields remain byte-identical after standard trailing-blank
padding. The production reader and checkpoint raw-history warmup both consume
this derived file; no production force-disable switch is introduced.
"""
import argparse
import hashlib
import json
from pathlib import Path


def sha(path):
    value=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):value.update(block)
    return value.hexdigest()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input-root',type=Path,required=True)
    parser.add_argument('--output-root',type=Path,required=True)
    args=parser.parse_args()
    name='MAL200KEN_R_20241990000_01D_30S_MO.rnx'
    source=args.input_root/'data'/name
    before=sha(source)
    args.output_root.mkdir(parents=True,exist_ok=False)
    (args.output_root/'data').mkdir()
    (args.output_root/'products').symlink_to((args.input_root/'products').resolve(),target_is_directory=True)
    target=args.output_root/'data'/name
    types={};counts={};system=None;header=True;records=0;cleared=0
    with source.open('rb') as reader,target.open('wb') as writer:
        for raw in reader:
            line=raw.rstrip(b'\r\n')
            if header:
                label=line[60:80].strip()
                if label==b'SYS / # / OBS TYPES':
                    if line[:1]!=b' ':
                        system=line[:1];counts[system]=int(line[3:6]);types[system]=[]
                    types[system].extend(line[7:60].split())
                if label==b'END OF HEADER':
                    if len(types.get(b'G',[]))!=counts.get(b'G'):raise RuntimeError('GPS type inventory incomplete')
                    header=False
                writer.write(raw);continue
            if line[:3]==b'   ' and line.strip():raise RuntimeError('Wrapped/event record unsupported; do not silently alter layout')
            if len(line)>=3 and line[:1]==b'G' and line[1:3].isdigit():
                names=types[b'G'];width=3+16*len(names)
                if len(line)>width:raise RuntimeError('Unexpected GPS record length')
                padded=line.ljust(width,b' ');modified=bytearray(padded)
                for index,observable in enumerate(names):
                    start=3+16*index
                    if observable.startswith(b'D'):
                        cleared+=bool(padded[start:start+16].strip())
                        modified[start:start+16]=b' '*16
                    elif modified[start:start+16]!=padded[start:start+16]:
                        raise RuntimeError('Non-Doppler observation changed')
                writer.write(bytes(modified)+b'\n');records+=1
            else:writer.write(raw)
    if header or not records or not cleared:raise RuntimeError('No complete GPS Doppler fixture produced')
    unchanged=sha(source)==before
    if not unchanged:raise RuntimeError('Original RINEX changed during fixture preparation')
    result={'original':str(source),'original_sha256':before,'derived':str(target),
        'derived_sha256':sha(target),'gps_records':records,'nonempty_doppler_fields_cleared':cleared,
        'gps_observable_types':[value.decode('ascii') for value in types[b'G']],
        'non_doppler_fields_identical':True,'original_unchanged':unchanged,
        'scope':'Derived no-Doppler QC input; original phase/code/LLI values are not injected or repaired'}
    (args.output_root/'fixture_manifest.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
