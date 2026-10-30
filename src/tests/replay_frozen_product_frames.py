#!/usr/bin/env python3
"""Convert frozen ProductFrame inputs without changing T/H/b or column order.

Generated physical IDs are synthetic labels, not authenticated arc identities.
This is a 02:00 exact-backend consistency experiment, NOT a 02:06 / 02:37
bridge replay, a covariance comparison, or an AR/product-quality validation.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--driver", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    converted = args.output / "converted_inputs"
    converted.mkdir(exist_ok=True)
    manifest = {
        "scope": "generic_real0200_inputs_exact_backend_consistency",
        "not_bridge_event_replay": True,
        "not_covariance_or_AR_search_validation": True,
        "physical_ids_are_synthetic": True,
        "preserves_original_column_order_and_sparse_payload": True,
        "query_work": "IMAGE_GENERATORS_ONLY existing capability",
        "driver_sha256": digest(args.driver),
        "inputs": [],
    }
    paths = sorted(args.inputs.glob("frame_*.txt"),
                   key=lambda p: int(p.stem.split("_")[-1]))
    if not paths:
        raise RuntimeError("no frozen product frame inputs")
    for number, source in enumerate(paths, 1):
        text = source.read_text(encoding="utf-8")
        header, payload = text.split("\n", 1)
        magic, dimension, work = header.split()
        if magic != "R51_PRODUCT_FRAME_V1":
            raise RuntimeError(f"unexpected frame format {source}")
        dimension = int(dimension)
        target = converted / source.name
        if target.exists():
            raise RuntimeError(f"refusing to overwrite {target}")
        prefix = [
            f"R51_BRIDGE_SUPPORT_V1 {dimension} {number}",
            json.dumps("2024-07-17 02:00:00") + " " +
            json.dumps(f"SYNTHETIC_FRAME_{number}") + " " +
            json.dumps("FROZEN_REAL_0200_ALGEBRA_ONLY"),
        ]
        prefix += [json.dumps(f"SYNTHETIC_COLUMN_{c}") for c in range(dimension)]
        prefix += ["0"]  # No real physical-arc availability evidence was saved.
        target.write_text("\n".join(prefix) + "\n" + payload, encoding="utf-8")
        manifest["inputs"].append({
            "source": str(source), "source_sha256": digest(source),
            "converted": str(target), "converted_sha256": digest(target),
            "original_work_enum": int(work), "dimension": dimension,
            "original_sparse_payload_sha256": hashlib.sha256(payload.encode()).hexdigest(),
        })
    (args.output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    command = ["/usr/bin/time", "-v", "-o", str(args.output / "time.txt"),
               str(args.driver), "--frozen-algebra", str(converted)]
    with (args.output / "results.jsonl").open("w", encoding="utf-8") as stdout, \
         (args.output / "stderr.log").open("w", encoding="utf-8") as stderr:
        result = subprocess.run(command, stdout=stdout, stderr=stderr, check=False)
    summary = {
        "returncode": result.returncode, "input_count": len(paths),
        "results_count": len((args.output / "results.jsonl").read_text().splitlines()),
        "command": command, "scope": manifest["scope"],
    }
    (args.output / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary))
    raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
