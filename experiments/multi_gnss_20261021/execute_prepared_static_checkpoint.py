"""Execute one still-PREPARED controlled plan with a separate pinned approval.

No configuration or physics are rewritten; solver command/provenance remain
identical to preparation. An existing/failed/finished run cannot be retried.
"""
from pathlib import Path
import json
import sys

import run_static_checkpoint_replay as replay
import run_stages as stages

path, approval = (Path(value).resolve(strict=True) for value in sys.argv[1:3])
plan = json.loads(path.read_text())
if plan.get("state") != "PREPARED_NOT_RUN":
    raise RuntimeError("Only an unused prepared controlled replay may be executed")
env = stages.environment("D")
env.update(plan["environment"])
try:
    replay.execute(plan, env, approval)
except Exception as error:
    plan.update(state="CONTROLLED_RESTORE_CONTRACT_FAILED", failure_reason=str(error),
                finished_utc=replay.stamp(), controlled_restore_contract_passed=False)
    stages.write_json(path, plan)
    raise
print(json.dumps({"state": plan["state"], "run_root": plan["run_root"]}), flush=True)
