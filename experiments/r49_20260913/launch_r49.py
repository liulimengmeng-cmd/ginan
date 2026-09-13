from pathlib import Path
import datetime,json,subprocess
a=Path(__file__).resolve().parent
spec=json.loads((a/'experiment.json').read_text())
frozen=json.loads((a/'frozen_r49.json').read_text())
with (a/'launch.log').open('xb') as log:
    process=subprocess.Popen(['bash',str(a/'run_r49_30m.sh')],stdin=subprocess.DEVNULL,
        stdout=log,stderr=subprocess.STDOUT,start_new_session=True,close_fds=True)
receipt={'runner_pid':process.pid,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'case':spec['case'],'source_commit':frozen['source_commit'],
    'frozen_root':frozen['frozen_root'],'binary_sha256':frozen['binary_sha256'],
    'observation_window':'2024-07-17 00:00:00 to 00:30:00','expected_epochs':61,
    'note':'Launch receipt only; actual PEA startup must be checked separately.'}
(a/'launch_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
