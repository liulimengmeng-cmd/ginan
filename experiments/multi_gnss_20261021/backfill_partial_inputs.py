#!/usr/bin/env python3
"""Fetch only same-name/day original-network files with actual epoch gaps.

An upstream file retaining the same gaps is recorded, never interpolated or
renamed to pretend completeness. Original files remain untouched.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip
import json
from pathlib import Path
import re
import shutil
import subprocess
from urllib.parse import urljoin

import download_daily_inputs as daily


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    args = parser.parse_args()
    root = args.root
    validation = json.loads((root / "validation_results.json").read_text())
    if validation.get("scanner_contract") != "rinex_fixed_width_space_padded_prn_and_undated_event_v2":
        raise RuntimeError("Require the corrected fixed-width scanner")
    if len(validation["existing_observations"]) != 180 or "summary" not in validation:
        raise RuntimeError("Wait for the full original 180-station audit")
    partial = [row for row in validation["existing_observations"] if not row["full_2880_epoch_grid"]]
    links = set(re.findall(r'href="([^"]+)"', (root / "source_indexes/observations.html").read_text()))
    audit = json.loads(args.audit.read_text())
    whitelist = {r["prn"] for r in audit["whitelist"]}
    manifest_path = root / "download_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    converter = next((root / "tools/RNXCMP_4.2.0").rglob("CRX2RNX.exe"))
    receipts = []
    jobs = []
    if len(partial) > 12:
        raise RuntimeError("Unexpectedly many partial files; do not silently expand data acquisition")
    for row in partial:
        name = Path(row["path"]).name.replace(".rnx", ".crx.gz")
        if name not in links:
            receipts.append({"original": row["path"], "status": "NOT_IN_BKG_DIRECTORY", "missing_epoch_count": len(row["missing_30s_epochs"])})
        else:
            jobs.append((row, urljoin(daily.OBS_URL, name), root / "observations" / name))
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures = {pool.submit(daily.get, url, path): (row, path) for row, url, path in jobs}
        for future in as_completed(futures):
            old, compressed = futures[future]
            download = future.result()
            download["purpose"] = "Same-station/day epoch-gap cross-check; never overwrite original"
            manifest["downloads"].append(download)
            daily.dump(manifest_path, manifest)
            print(json.dumps(download), flush=True)
            if download["status"] == "FAILED":
                receipts.append({"original": old["path"], "status": "DOWNLOAD_FAILED", "download": download})
                continue
            crx = compressed.with_suffix("")
            rnx = crx.with_suffix(".rnx")
            if not crx.exists():
                with gzip.open(compressed, "rb") as source, crx.open("xb") as output:
                    shutil.copyfileobj(source, output)
            conversion = {"converter": str(converter), "converter_sha256": daily.sha(converter)}
            if not rnx.exists():
                with crx.open("rb") as source, rnx.open("xb") as output:
                    process = subprocess.run([str(converter), "-"], stdin=source, stdout=output, stderr=subprocess.PIPE)
                conversion.update(exit_code=process.returncode, stderr=process.stderr.decode(errors="replace"))
                if process.returncode != 0:
                    receipts.append({"original": old["path"], "status": "CONVERSION_FAILED", "conversion": conversion})
                    continue
            result = daily.observation_audit(rnx, whitelist)
            receipt = {"station": rnx.name[:9], "original": old["path"], "downloaded_rinex": str(rnx),
                "old_epoch_count": old["epoch_count"], "new_epoch_count": result["epoch_count"],
                "old_missing_seconds": old["missing_30s_epochs"], "new_missing_seconds": result["missing_30s_epochs"],
                "normalized_observation_body_sha_equal": old["normalized_observation_body_sha256"] == result["normalized_observation_body_sha256"],
                "status": "NEW_COPY_COMPLETE" if result["full_2880_epoch_grid"] else "UPSTREAM_SAME_DAY_COPY_STILL_PARTIAL_NO_FABRICATED_FILL",
                "download": download, "conversion": conversion, "audit": result}
            receipts.append(receipt)
            daily.dump(root / "partial_backfill_audit.json", receipts)
    daily.dump(root / "partial_backfill_audit.json", receipts)
    print(json.dumps([{k: row.get(k) for k in ("station", "status", "old_epoch_count", "new_epoch_count", "new_missing_seconds", "normalized_observation_body_sha_equal")} for row in receipts], indent=2), flush=True)


if __name__ == "__main__":
    main()
