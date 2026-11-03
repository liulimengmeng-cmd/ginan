"""Generate only the companion-declared private DTO mirror; never compile/run."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def render(root, companion_sha):
    companion=Path(root)/'src/tests/run_gec_ar_runtime_checkpoint.py'
    if digest(companion)!=companion_sha:
        raise RuntimeError('Exact frozen companion SHA required')
    text=companion.read_text()
    module=ast.parse(text)
    nodes=[n for n in module.body if isinstance(n,ast.FunctionDef) and n.name=='declaration']
    if len(nodes)!=1:
        raise RuntimeError('Exactly one original declaration extraction function required')
    function=ast.get_source_segment(text,nodes[0])
    source_path=Path(root)/'src/cpp/pea/ppp_ambres_checkpoint.inc'
    source=source_path.read_text().split('/** Multi-GNSS FLOAT checkpoints')[0]
    scope={'re':re,'source':source}
    exec(compile(function,str(companion)+':declaration-only','exec'),scope)
    names=re.findall(r'^struct (\w+)\s*\{',source,re.MULTILINE)
    if not names or len(names)!=len(set(names)):
        raise RuntimeError('Private DTO names missing or duplicated')
    mirror='#pragma once\nnamespace gec_test_ambres_wire {\n'+\
        '\n\n'.join(scope['declaration'](name) for name in names)+'\n}\n'
    proof=dict(schema='GEC_PRIVATE_DTO_DECLARATIONS_ONLY_V1',
        source_path=str(source_path),source_sha256=digest(source_path),
        companion_path=str(companion),companion_sha256=companion_sha,
        extraction_function_sha256=hashlib.sha256(function.encode()).hexdigest(),
        declaration_names=names,header_sha256=hashlib.sha256(mirror.encode()).hexdigest(),
        validators_or_exporters_copied=False,no_compiler_or_solver_run=True)
    return mirror,proof


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-root',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--companion-sha256',required=True)
    args=p.parse_args()
    if args.output_dir.exists():
        raise RuntimeError('Fresh schema staging directory required; preserve prior evidence')
    mirror,proof=render(args.source_root,args.companion_sha256)
    args.output_dir.mkdir(parents=True)
    header=args.output_dir/'gec_private_ambres_wire_schema.hpp'
    with header.open('x',encoding='utf-8',newline='\n') as f:f.write(mirror)
    with (args.output_dir/'PROVENANCE.json').open('x',encoding='utf-8',newline='\n') as f:
        json.dump(proof,f,indent=2);f.write('\n')
    print(json.dumps(dict(generated=True,header=str(header),sha256=digest(header),
        declarations=len(proof['declaration_names']),no_compiler_or_solver_run=True),indent=2))


if __name__=='__main__':
    main()
