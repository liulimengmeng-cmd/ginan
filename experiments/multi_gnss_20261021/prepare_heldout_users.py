#!/usr/bin/env python3
"""Prepare paired independent STATIC users; never start a solver in this tool.

F remains blocked until real G/E/C products AND independently validated
multi-system GONSIS/clock/CSV65/full-Q consumer support exist. Config generation
must never stand in for those results. Original user configurations are sealed.
"""
from __future__ import annotations

import argparse
import copy
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

import yaml

import run_stages

CLIENT_MANIFEST = "artifacts/user_pro_product_model_20261009/experiments_static_isolate_20261009T125203Z/manifest.json"
SIGNALS = {"gps": ["L1C", "L2W"], "gal": ["L1C", "L5Q"], "bds": ["L1P", "L5P"]}
SELECTED = ["REYK", "PTGG"]


def sha(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def receipt(path):
    return {"path": str(path), "size_bytes": path.stat().st_size, "sha256": sha(path)}


def native(path):
    value = str(path).replace("\\", "/")
    if re.match("^[A-Za-z]:/", value):
        return "/mnt/" + value[0].lower() + value[2:]
    return value


def basename(path):
    return str(path).replace("\\", "/").rsplit("/", 1)[-1]


def host_path(value, work):
    """Resolve sealed /mnt/c/... receipts on Windows or Linux, without guessing."""
    if re.match(r"^/mnt/[a-z]/", value) and re.match(r"^[A-Za-z]:", str(work)):
        return Path(value[5].upper() + ":/" + value[7:])
    return Path(value)


def reference_coordinates(path, stations):
    block, selected = "", {}
    with path.open("rb") as stream:
        for number, raw in enumerate(stream, 1):
            line = raw.decode("ascii", "replace").rstrip("\r\n")
            if line.startswith("+"):
                block = line.strip()
            elif line.startswith("-"):
                block = ""
            if block != "+SOLUTION/ESTIMATE" or not line.startswith(" "):
                continue
            fields = line.split()
            if len(fields) < 10 or fields[1] not in ("STAX", "STAY", "STAZ") or fields[2] not in stations:
                continue
            if fields[5] != "24:199:43200":
                continue
            selected.setdefault(fields[2], {})[fields[1]] = {
                "estimate_m": float(fields[8]), "sigma_m": float(fields[9]),
                "point": fields[3], "solution": fields[4], "reference_epoch": fields[5],
                "source_line": number, "source_line_sha256": hashlib.sha256(raw).hexdigest()}
    for station in stations:
        if set(selected.get(station, {})) != {"STAX", "STAY", "STAZ"}:
            raise RuntimeError("Missing independent SINEX XYZ for " + station)
        if len({(row["point"], row["solution"]) for row in selected[station].values()}) != 1:
            raise RuntimeError("Mixed SINEX point/solution XYZ")
    return {"source": receipt(path), "type": "IGS weekly combined SINEX SOLUTION/ESTIMATE, not APRIORI",
            "frame": "IGS20", "epoch": "2024-07-17 12:00:00", "coordinates": selected,
            "used_as_user_coordinate_constraint": False,
            "warning": "External coordinate comparison is independent of this 12-station service; it is not independent orbit/clock truth."}


def require_heldout(station, readiness, validation):
    if station in readiness["selected_12_station_candidates"]:
        raise RuntimeError(station + " is in the service network")
    rows = [row for row in validation["existing_observations"] if basename(row["path"])[:4] == station]
    if len(rows) != 1:
        raise RuntimeError(station + " is not uniquely in the original 180-station input inventory")
    row = rows[0]
    headers = [h for h in readiness["station_headers"] if h["filename"][:4] == station]
    if len(headers) != 1 or not all(headers[0]["header_target_support"].values()) or not all(headers[0]["receiver_antex_target_frequency_support"].values()):
        raise RuntimeError(station + " lacks exact G/E/C signals or receiver ANTEX frequencies")
    if not row["full_2880_epoch_grid"] or row["parse_error_count"] or any(row["target_nonzero_pair_epoch_counts"][s] != 2880 for s in "GEC"):
        raise RuntimeError(station + " lacks full-day raw target-pair observations")
    return row


def configuration(template, station, mode, observation, output_root, products, covariance, readiness):
    cfg = copy.deepcopy(template)
    case = "multi_gnss_F_" + station + "_static_" + mode
    cfg["inputs"]["inputs_root"] = "/mnt/d/GINAN_R20/inputData/products/"
    cfg["inputs"]["gnss_observations"] = {"gnss_observations_root": "/mnt/d/GINAN_R20/inputData/data/",
                                          "rnx_inputs": [basename(observation["path"])]}
    # Truth coordinates are a scoring input only. Never add the coordinate
    # SINEX to the solver input or fix the user station to its truth.
    if any("20241960000_07D_07D_CRD" in value for value in cfg["inputs"].get("snx_files", [])):
        raise RuntimeError("Unexpected truth-coordinate input in the inherited user model")
    cfg["outputs"]["outputs_root"] = output_root + "/outputs/" + case
    cfg["outputs"]["metadata"]["config_description"] = case
    for output_kind in ("bias_sinex", "clocks", "sp3", "orbex", "sinex"):
        if output_kind in cfg["outputs"]:
            cfg["outputs"][output_kind]["output"] = False
    gnss = cfg["processing_options"]["gnss_general"]
    gnss.update(use_tgd_bias=False, minimise_ionosphere_offsets=False)
    gnss["zhang_full_rank"] = {"enable": False, "bds3_only": True,
        "bds3_satellite_whitelist": [row["prn"] for row in readiness["whitelist"]]}
    gnss["sys_options"] = {}
    for system in ("gps", "gal", "bds", "glo", "qzs", "sbs", "irn", "leo"):
        options = {"process": system in SIGNALS, "ambiguity_resolution": system in SIGNALS}
        if system in SIGNALS:
            options.update(code_priorities=SIGNALS[system], receiver_amb_pivot=True,
                           network_amb_pivot=False, reject_eclipse=False,
                           used_nav_type=run_stages.NAV_TYPES[system])
        gnss["sys_options"][system] = options
    z = gnss["zhang_pppar"]
    z.update(user_adapter=True, output_products=False, deterministic_checkpoint=False,
             checkpoint_runtime_id=case, product_filename=products, product_covariance_filename=covariance,
             product_solution="PRODUCT_FIXED", user_use_full_product_covariance=True,
             sys_options={s: {"baseline_observables": pair} for s, pair in SIGNALS.items()})
    cfg["processing_options"]["ambiguity_resolution"].update(
        mode="OFF" if mode == "FLOAT" else "LAMBDA_ALT", solution_ratio_threshold=2,
        elevation_mask=15, once_per_epoch=True, fix_and_hold=False, lambda_set_size=2)
    cfg["processing_options"]["epoch_control"].update(start_epoch="2024-07-17 02:00:00",
        end_epoch="2024-07-17 02:59:30", epoch_interval=30)
    cfg["receiver_options"] = {"global": cfg["receiver_options"]["global"]}
    cfg["receiver_options"]["global"]["elevation_mask"] = 7.5
    cfg["estimation_parameters"]["receivers"] = {"global": cfg["estimation_parameters"]["receivers"]["global"]}
    cfg["estimation_parameters"]["receivers"]["global"]["pos"].update(estimated=[True], process_noise=[0])
    satellite_options = run_stages.overlay("E", Path("/unused"), readiness)["satellite_options"]
    satellite_options["global"] = cfg["satellite_options"]["global"]
    if "GPS" in cfg["satellite_options"]:
        satellite_options["GPS"] = cfg["satellite_options"]["GPS"]
    cfg["satellite_options"] = satellite_options
    return cfg


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--readiness", type=Path, required=True)
    parser.add_argument("--validation", type=Path, required=True)
    parser.add_argument("--coordinate-sinex", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--native-user-run-root", required=True)
    parser.add_argument("--service-run-root", required=True)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    (out / "configs").mkdir()
    readiness = json.loads(args.readiness.read_text())
    validation = json.loads(args.validation.read_text())
    coordinates = reference_coordinates(args.coordinate_sinex, SELECTED)
    source_manifest_path = args.work / CLIENT_MANIFEST
    source = json.loads(source_manifest_path.read_text())
    templates = {}
    source_cases = {}
    for mode, strategy in (("FLOAT", "FLOAT"), ("AR", "GONSIS_UDUC")):
        selected = next(case for case in source["cases"] if case["station"] == "DYNG" and
                        case["provider"] == "internal" and case["position_mode"] == "static" and case["strategy"] == strategy)
        path = host_path(selected["config"]["path"], args.work)
        if sha(path) != selected["config"]["sha256"]:
            raise RuntimeError("Sealed GONSIS source configuration changed")
        templates[mode] = yaml.safe_load(path.read_text())
        source_cases[mode] = dict(config=receipt(path), environment=selected["environment"])
    configs = []
    # Consume the service's immutable per-epoch bundles directly. An aggregate
    # CSV is deliberately not appended by the atomic writer and is not authority.
    exported_product = args.service_run_root.rstrip("/") + "/outputs/zhang_internal_products.csv"
    exported_covariance = args.service_run_root.rstrip("/") + "/outputs/zhang_internal_product_covariance.csv"
    for station in SELECTED:
        observation = require_heldout(station, readiness, validation)
        for mode in ("FLOAT", "AR"):
            cfg = configuration(templates[mode], station, mode, observation, args.native_user_run_root,
                                exported_product, exported_covariance, readiness)
            path = out / "configs" / (station + "_static_" + mode + ".yaml")
            path.write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")
            env = copy.deepcopy(source_cases[mode]["environment"])
            env.update(ZHANG_R51_RATIO_ONLY="1", ZHANG_USER_FLOAT_TRANSPORT="CONTRACT",
                       ZHANG_USER_UNKNOWN_CHART="ISOLATE", ZHANG_USER_GONSIS_TARGETED_PAR="1",
                       ZHANG_USER_AR_FEEDBACK="EPOCH_ONLY", ZHANG_USER_PRODUCT_ERROR_POLICY="WHITE",
                       ZHANG_USER_PRODUCT_ERROR_UNSUPPORTED="FAIL", ZHANG_USER_SYSTEM_CLOCKS="1")
            configs.append({"station": station, "mode": mode, "config": receipt(path),
                "native_config": native(path), "environment": env,
                "observation": {k: observation[k] for k in ("path", "sha256", "size_bytes", "full_2880_epoch_grid", "target_nonzero_pair_epoch_counts")},
                "truth_coordinate_input": False})
    # Paired runs differ only by AR mode/output identity. Enforce same data,
    # stochastic model, gates and product errors before a future real launch.
    for station in SELECTED:
        float_cfg = json.loads((out / "configs" / (station + "_static_FLOAT.yaml")).read_text())
        ar_cfg = json.loads((out / "configs" / (station + "_static_AR.yaml")).read_text())
        for cfg in (float_cfg, ar_cfg):
            cfg["outputs"] = {}
            cfg["processing_options"]["ambiguity_resolution"]["mode"] = "NORMALIZED_FOR_PAIR_CHECK"
            cfg["processing_options"]["gnss_general"]["zhang_pppar"]["checkpoint_runtime_id"] = "PAIR"
        if float_cfg != ar_cfg:
            raise RuntimeError("Unexpected FLOAT/AR physical model mismatch")
    excluded = {station: "Not in original 180-station service inventory" for station in ("DYNG", "NICO")}
    manifest = {"status": "BLOCKED_NOT_RUN", "prepared_utc": datetime.now(timezone.utc).isoformat(),
        "solver_launched": False, "scientific_acceptance": False, "selected_stations": SELECTED,
        "service_stations": readiness["selected_12_station_candidates"], "preferred_station_exclusions": excluded,
        "source_manifest": receipt(source_manifest_path), "source_cases": source_cases,
        "reference_coordinates": coordinates, "readiness": receipt(args.readiness), "validation": receipt(args.validation),
        "cases": configs, "service_run_root": args.service_run_root,
        "native_user_run_root": args.native_user_run_root, "expected_product_epochs": 120,
        "publication_contract": {"minimum_columns": 65, "preserve_all_original_flags": True,
            "authority": "Immutable R51 v2 epoch bundle with paired covariance and certificate",
            "same_segment_continuous": "Use published proof; never infer from smooth phase values",
            "requires_user_ambiguity_reset": "Obey actual contract; a new segment is not silently treated as physical continuity",
            "known_transport": "CONTRACT", "unknown_transport": "ISOLATE", "never_fill_cross_covariance_with_zero": True},
        "blocking_gates": ["Real committed G/E/C service products for 02:00:00-02:59:30 are not available yet",
            "Multi-system GONSIS WL/L1/PAR consumer is not scientifically validated by the older GPS-only evidence",
            "Independent user effective clocks and code/clock datum compatibility must be checked in the actual user binary",
            "Elementwise full Q reconstruction, nonzero G/E-G/C-E/C cross blocks, PSD and shared epoch-factor identity must pass",
            "CSV65 continuity contract and E11 slip locality must pass without fabricating any flags"],
        "not_a_result": "Prepared configs do not establish accuracy, convergence or correct-fix rate"}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": manifest["status"], "selected_stations": SELECTED,
                      "configs": len(configs), "solver_launched": False}, indent=2))


if __name__ == "__main__":
    main()
