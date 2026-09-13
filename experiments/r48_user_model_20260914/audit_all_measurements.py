from pathlib import Path
w=Path('/mnt/c/Users/rx/Documents/GINAN/r48_user_model_work_20260914')
src=(w/'independent_posterior_audit.py').read_text()
start=src.index('  d,label=read(files[-1])');end=src.index("(w/'independent_posterior_audit.json')")
block=src[start:end].replace('read(files[-1])','read(file)')
src=src[:start]+'  for file in files:\n'+''.join(' '+line+'\n' for line in block.splitlines())+src[end:].replace('independent_posterior_audit.json','independent_all_measurements.json')
exec(compile(src,'independent_all_measurements','exec'))
