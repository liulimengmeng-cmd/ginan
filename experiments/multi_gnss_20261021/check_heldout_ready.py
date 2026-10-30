#!/usr/bin/env python3
"""Fail-closed F preflight. This command never launches a scientific process."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from decimal import Decimal

from prepare_heldout_users import sha


REQUIRED = ("build_passed", "mandatory_tests_passed", "multi_system_user_search_validated",
            "independent_user_clock_model_verified", "precise_clock_datum_review_passed",
            "csv65_continuity_contract_verified", "full_joint_covariance_verified",
            "cross_system_shared_epoch_factor_verified", "e11_slip_isolation_verified",
            "atomic_v2_bundle_consumer_verified")


def audit_atomic_bundles(base):
    """Independent structural/hash preflight, not an integer-proof verifier."""
    from run_stages import epoch_seconds
    expected = {epoch_seconds("02:00:00") + 30 * i for i in range(120)}
    records = {}
    if not base.is_dir():
        raise RuntimeError("NO_REAL_MULTI_SYSTEM_ATOMIC_SERVICE_BUNDLES")
    for directory in sorted(base.iterdir()):
        if not directory.is_dir() or directory.name.endswith(".pending"):
            continue
        commit = directory / "COMMIT.json"
        if not commit.is_file():
            continue
        metadata = json.loads(commit.read_text())
        epoch_text = metadata.get("epoch", "")
        if not epoch_text.startswith("2024-07-17 "):
            continue
        epoch = epoch_seconds(epoch_text.split(" ", 1)[1])
        if epoch not in expected:
            continue
        if epoch in records:
            raise RuntimeError("AMBIGUOUS_SERVICE_EPOCH")
        if metadata.get("schema_version") != 2 or set(metadata.get("systems", {})) != {"GPS", "GAL", "BDS"}:
            raise RuntimeError("MULTI_SYSTEM_V2_AUTHORITY_MISSING")
        names = {"COMMIT.json", "products.csv", "covariance.csv", "certificate.txt"}
        sums = {}
        for line in (directory / "SHA256SUMS").read_text().splitlines():
            digest, name = line.split(None, 1)
            if name not in names or name in sums:
                raise RuntimeError("SERVICE_HASH_MANIFEST_INVALID")
            sums[name] = digest
        if set(sums) != names:
            raise RuntimeError("SERVICE_HASH_MANIFEST_INCOMPLETE")
        for name, digest in sums.items():
            path = directory / name
            if not path.is_file() or sha(path) != digest:
                raise RuntimeError("SERVICE_BUNDLE_HASH_MISMATCH:" + name)
        if (directory / "products.csv").stat().st_size != metadata["products_bytes"] or \
                (directory / "covariance.csv").stat().st_size != metadata["covariance_bytes"]:
            raise RuntimeError("SERVICE_BUNDLE_SIZE_MISMATCH")
        with (directory / "products.csv").open(newline="") as stream:
            reader = csv.DictReader(stream)
            required = {"same_segment_continuous", "requires_user_ambiguity_reset",
                        "public_phase_segment_id", "joint_covariance_complete"}
            if len(reader.fieldnames or []) != 65 or not required <= set(reader.fieldnames or []):
                raise RuntimeError("PRODUCT_LACKS_PUBLIC65_CONTINUITY_CONTRACT")
            systems = set()
            for row in reader:
                if Decimal(row["gpst_seconds"]) != epoch:
                    raise RuntimeError("SERVICE_CSV_COMMIT_EPOCH_MISMATCH")
                if row["solution"] == "PRODUCT_FIXED":
                    systems.add(row["satellite"][0])
            if systems != {"G", "E", "C"}:
                raise RuntimeError("SERVICE_EPOCH_GEC_PRODUCT_COVERAGE_MISSING")
        records[epoch] = sha(commit) + ":" + sha(directory / "SHA256SUMS")
    if set(records) != expected:
        raise RuntimeError("SERVICE_ATOMIC_120_EPOCH_COVERAGE_MISSING:" + str(len(records)))
    return hashlib.sha256(json.dumps(records, sort_keys=True).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--acceptance", type=Path)
    parser.add_argument("--binary", type=Path)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    blocks = []
    acceptance = json.loads(args.acceptance.read_text()) if args.acceptance else {}
    for field in REQUIRED:
        if acceptance.get(field) is not True:
            blocks.append("MISSING_ACTUAL_ACCEPTANCE:" + field)
    if args.binary is None or not args.binary.is_file():
        blocks.append("NO_NEW_MULTI_SYSTEM_USER_BINARY")
    elif acceptance.get("binary_sha256") != sha(args.binary):
        blocks.append("USER_BINARY_ACCEPTANCE_SHA256_MISMATCH")
    for case in manifest["cases"]:
        path = Path(case["native_config"])
        if not path.is_file() or sha(path) != case["config"]["sha256"]:
            blocks.append("PREPARED_CONFIG_CHANGED_OR_MISSING:" + case["station"] + ":" + case["mode"])
    first = manifest["cases"][0]
    config = json.loads(Path(first["native_config"]).read_text())
    z = config["processing_options"]["gnss_general"]["zhang_pppar"]
    try:
        fingerprint = audit_atomic_bundles(Path(z["product_filename"] + ".epochs"))
        if acceptance.get("atomic_bundle_set_sha256") != fingerprint:
            blocks.append("SERVICE_ATOMIC_BUNDLE_SET_NOT_PINNED")
    except (OSError, ValueError, KeyError, RuntimeError) as error:
        blocks.append(str(error))
    print(json.dumps({"status": "READY_FOR_SEPARATE_EXPLICIT_LAUNCH" if not blocks else "BLOCKED_NOT_RUN",
                      "solver_launched": False, "scientific_acceptance": False, "blocks": blocks}, indent=2))
    raise SystemExit(1 if blocks else 0)


if __name__ == "__main__":
    main()
