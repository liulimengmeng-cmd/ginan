#!/usr/bin/env python3
"""Freeze and evaluate the Experiment 0 2024 float-covariance pivot.

The original fixed-STEC prerequisite failed on both registered days.  This
module therefore implements a deliberately narrower, leakage-resistant pivot:

* model and validation targets both come from FLOAT forward-filter sidecars;
* the response is the datum-safe satellite-minus-reference STEC contrast;
* the mean is a fixed thin-shell VTEC plane;
* ``none``, ``diagonal`` and ``full`` use the same rows and mean model;
* quiet-day model-station leave-one-site-out residuals calibrate one variance
  inflation factor per arm before any held-out output is read;
* held-out FLOAT covariance is included only as validation-target uncertainty.

This is statistical-pipeline evidence.  It cannot establish fixed STEC,
integer correctness, an independent geophysical gradient truth, or the
original PPP-AR wrong-fix hypotheses.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from decimal import Decimal
from pathlib import Path
from typing import Iterable, Iterator, Sequence, TextIO

import numpy as np

from validate_stec_satellite_difference_covariance import (
    DifferenceEpoch,
    iter_difference_epochs,
)


FREEZE_SCHEMA = "GINAN_EXPERIMENT0_2024_FLOAT_SPATIAL_FREEZE_V1"
REPORT_SCHEMA = "GINAN_EXPERIMENT0_2024_FLOAT_SPATIAL_VALIDATION_V1"
METHODS = ("none", "diagonal", "full")
GPS_EPOCH = datetime(1980, 1, 6, tzinfo=timezone.utc)
EARTH_RADIUS_M = 6_371_000.0
SHELL_HEIGHT_M = 506_700.0
MSLM_ALPHA = 0.9782
ORIGIN_LAT_RAD = math.radians(-33.0)
ORIGIN_LON_RAD = math.radians(145.0)
MIN_ELEVATION_RAD = math.radians(15.0)
WGS84_A = 6_378_137.0
WGS84_F = 1 / 298.257223563
WGS84_E2 = WGS84_F * (2 - WGS84_F)
REGULARISATION_ABSOLUTE = 1e-12
REGULARISATION_RELATIVE = 1e-10
BOOTSTRAP_REPLICATES = 1000
BOOTSTRAP_SEED = 20240508


class SpatialInputError(ValueError):
    """Raised when frozen inputs cannot be interpreted without ambiguity."""


@dataclass(frozen=True)
class RayGeometry:
    site: str
    satellite: str
    state_number: int
    elevation_rad: float
    ipp_lat_rad: float
    ipp_lon_rad: float
    mapping_factor: float
    feature: np.ndarray


@dataclass
class SpatialEpoch:
    gps_week: int
    gps_tow: Decimal
    status: str
    posterior_stage: str
    sites: np.ndarray
    target_satellites: list[str]
    reference_satellites: list[str]
    design: np.ndarray
    values: np.ndarray
    covariance: np.ndarray
    declared_difference_count: int
    dropped_geometry_count: int

    @property
    def key(self) -> tuple[int, Decimal]:
        return self.gps_week, self.gps_tow


@dataclass(frozen=True)
class FitResult:
    method: str
    beta: np.ndarray
    beta_covariance: np.ndarray
    rank: int
    training_count: int
    residual_variance: float
    covariance_rank: int


@dataclass
class CalibrationRecord:
    gps_week: int
    gps_tow: Decimal
    site: str
    residual: float
    prediction_variance: float
    target_variance: float

    @property
    def block(self) -> str:
        return f"{self.gps_week}:{int(self.gps_tow // Decimal(3600))}"


@dataclass
class PredictionRow:
    gps_week: int
    gps_tow: Decimal
    site: str
    target_satellite: str
    reference_satellite: str
    observed_tecu: float
    methods: dict[str, dict[str, float]]

    @property
    def block(self) -> str:
        return f"{self.gps_week}:{int(self.gps_tow // Decimal(3600))}"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _atomic_write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as stream:
            json.dump(payload, stream, indent=2, allow_nan=False)
            stream.write("\n")
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def _epoch_key(week: int, tow_text: str) -> tuple[int, Decimal]:
    tow = Decimal(tow_text)
    if week < 0 or not tow.is_finite() or tow < 0 or tow >= Decimal(604800):
        raise SpatialInputError(f"invalid GPS epoch {week}:{tow_text}")
    return week, tow


def _ecef_to_geodetic(ecef: np.ndarray) -> tuple[float, float, float]:
    x, y, z = (float(value) for value in ecef)
    longitude = math.atan2(y, x)
    horizontal = math.hypot(x, y)
    latitude = math.atan2(z, horizontal * (1 - WGS84_E2))
    height = 0.0
    for _ in range(12):
        sin_latitude = math.sin(latitude)
        prime_vertical = WGS84_A / math.sqrt(
            1 - WGS84_E2 * sin_latitude * sin_latitude
        )
        height = horizontal / max(math.cos(latitude), 1e-15) - prime_vertical
        next_latitude = math.atan2(
            z,
            horizontal
            * (1 - WGS84_E2 * prime_vertical / (prime_vertical + height)),
        )
        if abs(next_latitude - latitude) < 1e-14:
            latitude = next_latitude
            break
        latitude = next_latitude
    return latitude, longitude, height


def _azimuth_elevation(
    receiver_ecef: np.ndarray, satellite_ecef: np.ndarray
) -> tuple[float, float, float, float]:
    latitude, longitude, _ = _ecef_to_geodetic(receiver_ecef)
    delta = satellite_ecef - receiver_ecef
    sin_latitude = math.sin(latitude)
    cos_latitude = math.cos(latitude)
    sin_longitude = math.sin(longitude)
    cos_longitude = math.cos(longitude)
    east = -sin_longitude * delta[0] + cos_longitude * delta[1]
    north = (
        -sin_latitude * cos_longitude * delta[0]
        - sin_latitude * sin_longitude * delta[1]
        + cos_latitude * delta[2]
    )
    up = (
        cos_latitude * cos_longitude * delta[0]
        + cos_latitude * sin_longitude * delta[1]
        + sin_latitude * delta[2]
    )
    azimuth = math.atan2(east, north) % (2 * math.pi)
    elevation = math.atan2(up, math.hypot(east, north))
    return latitude, longitude, azimuth, elevation


def ray_geometry(
    site: str,
    satellite: str,
    state_number: int,
    receiver_ecef: np.ndarray,
    satellite_ecef: np.ndarray,
) -> RayGeometry:
    latitude, longitude, azimuth, elevation = _azimuth_elevation(
        receiver_ecef, satellite_ecef
    )
    shell_radius = EARTH_RADIUS_M + SHELL_HEIGHT_M
    rp = EARTH_RADIUS_M / shell_radius * math.cos(elevation)
    rp = min(1.0, max(-1.0, rp))
    central_angle = math.pi / 2 - elevation - math.asin(rp)
    sin_latitude = math.sin(latitude)
    cos_latitude = math.cos(latitude)
    sin_central = math.sin(central_angle)
    cos_central = math.cos(central_angle)
    ipp_latitude = math.asin(
        sin_latitude * cos_central
        + cos_latitude * sin_central * math.cos(azimuth)
    )
    longitude_increment = math.atan2(
        math.sin(azimuth) * sin_central * cos_latitude,
        cos_central - sin_latitude * math.sin(ipp_latitude),
    )
    ipp_longitude = (longitude + longitude_increment + math.pi) % (2 * math.pi) - math.pi

    map_argument = (
        EARTH_RADIUS_M
        / shell_radius
        * math.sin(MSLM_ALPHA * (math.pi / 2 - elevation))
    )
    mapping_factor = 1 / math.sqrt(max(1e-15, 1 - map_argument * map_argument))

    delta_longitude = (ipp_longitude - ORIGIN_LON_RAD + math.pi) % (2 * math.pi) - math.pi
    east_1000km = shell_radius * math.cos(ORIGIN_LAT_RAD) * delta_longitude / 1e6
    north_1000km = shell_radius * (ipp_latitude - ORIGIN_LAT_RAD) / 1e6
    feature = mapping_factor * np.array([1.0, east_1000km, north_1000km])
    return RayGeometry(
        site=site,
        satellite=satellite,
        state_number=state_number,
        elevation_rad=elevation,
        ipp_lat_rad=ipp_latitude,
        ipp_lon_rad=ipp_longitude,
        mapping_factor=mapping_factor,
        feature=feature,
    )


def iter_geometry_epochs(
    stream: TextIO,
) -> Iterator[tuple[tuple[int, Decimal], dict[tuple[str, str, int], RayGeometry]]]:
    current_key: tuple[int, Decimal] | None = None
    current: dict[tuple[str, str, int], RayGeometry] = {}
    for line_number, raw_line in enumerate(stream, start=1):
        if not raw_line.strip() or raw_line.startswith("#"):
            continue
        row = next(csv.reader([raw_line]))
        if len(row) != 16 or row[0] != "IONO_MEA":
            raise SpatialInputError(
                f"line {line_number}: expected 16-field IONO_MEA record"
            )
        key = _epoch_key(int(row[1]), row[2])
        if current_key is not None and key < current_key:
            raise SpatialInputError("STEC geometry epochs are not monotonic")
        if current_key is not None and key != current_key:
            yield current_key, current
            current = {}
        current_key = key
        site = row[3]
        satellite = row[4]
        state_number = int(row[7])
        geometry_key = (site, satellite, state_number)
        if geometry_key in current:
            raise SpatialInputError(
                f"duplicate STEC geometry row at {key}: {geometry_key}"
            )
        satellite_ecef = np.array([float(row[9]), float(row[10]), float(row[11])])
        receiver_ecef = np.array([float(row[13]), float(row[14]), float(row[15])])
        geometry = ray_geometry(
            site, satellite, state_number, receiver_ecef, satellite_ecef
        )
        if geometry.elevation_rad >= MIN_ELEVATION_RAD:
            current[geometry_key] = geometry
    if current_key is not None:
        yield current_key, current


def _covariance_matrix(epoch: DifferenceEpoch) -> np.ndarray:
    count = epoch.difference_state_count
    expected = count * (count + 1) // 2
    if len(epoch.covariance) != expected:
        raise SpatialInputError(
            f"incomplete SD covariance at {epoch.gps_week}:{epoch.gps_tow_text}"
        )
    covariance = np.full((count, count), np.nan)
    for (row, column), value in epoch.covariance.items():
        covariance[row, column] = value
        covariance[column, row] = value
    if not np.isfinite(covariance).all():
        raise SpatialInputError("nonfinite SD covariance")
    return covariance


def iter_spatial_epochs(stec_path: Path, sd_path: Path) -> Iterator[SpatialEpoch]:
    with stec_path.open("r", encoding="utf-8", errors="strict", newline="") as stec_stream:
        geometry_iterator = iter_geometry_epochs(stec_stream)
        geometry_item = next(geometry_iterator, None)
        with sd_path.open("r", encoding="utf-8", errors="strict", newline="") as sd_stream:
            for epoch in iter_difference_epochs(sd_stream):
                key = _epoch_key(epoch.gps_week, epoch.gps_tow_text)
                while geometry_item is not None and geometry_item[0] < key:
                    raise SpatialInputError(
                        f"geometry epoch {geometry_item[0]} has no matching SD epoch"
                    )
                geometry = {}
                if geometry_item is not None and geometry_item[0] == key:
                    geometry = geometry_item[1]
                    geometry_item = next(geometry_iterator, None)

                count = epoch.difference_state_count
                if sorted(epoch.states) != list(range(count)):
                    raise SpatialInputError(f"noncontiguous SD state catalogue at {key}")
                covariance = _covariance_matrix(epoch)
                kept: list[int] = []
                design_rows: list[np.ndarray] = []
                sites: list[str] = []
                targets: list[str] = []
                references: list[str] = []
                values: list[float] = []
                for local_index in range(count):
                    state = epoch.states[local_index]
                    target = geometry.get(
                        (state.site, state.target_satellite, state.state_number)
                    )
                    reference = geometry.get(
                        (state.site, state.reference_satellite, state.state_number)
                    )
                    if target is None or reference is None:
                        continue
                    kept.append(local_index)
                    design_rows.append(target.feature - reference.feature)
                    sites.append(state.site)
                    targets.append(state.target_satellite)
                    references.append(state.reference_satellite)
                    values.append(state.estimate_tecu)

                index = np.array(kept, dtype=int)
                subset = covariance[np.ix_(index, index)] if kept else np.zeros((0, 0))
                yield SpatialEpoch(
                    gps_week=epoch.gps_week,
                    gps_tow=key[1],
                    status=epoch.status,
                    posterior_stage=epoch.posterior_stage,
                    sites=np.array(sites, dtype=object),
                    target_satellites=targets,
                    reference_satellites=references,
                    design=np.vstack(design_rows) if design_rows else np.zeros((0, 3)),
                    values=np.array(values, dtype=float),
                    covariance=subset,
                    declared_difference_count=count,
                    dropped_geometry_count=count - len(kept),
                )
        if geometry_item is not None:
            raise SpatialInputError(
                f"trailing geometry epoch {geometry_item[0]} has no SD epoch"
            )


def _symmetric_pseudoinverse(matrix: np.ndarray) -> tuple[np.ndarray, int]:
    symmetric = (matrix + matrix.T) / 2
    eigenvalues, eigenvectors = np.linalg.eigh(symmetric)
    scale = max(1.0, float(np.max(np.abs(symmetric))))
    tolerance = REGULARISATION_ABSOLUTE + REGULARISATION_RELATIVE * scale
    if eigenvalues.size and float(eigenvalues[0]) < -tolerance:
        raise SpatialInputError(
            f"matrix is not PSD: min={eigenvalues[0]}, tolerance={tolerance}"
        )
    positive = eigenvalues > tolerance
    inverse = np.zeros_like(eigenvalues)
    inverse[positive] = 1 / eigenvalues[positive]
    return (eigenvectors * inverse) @ eigenvectors.T, int(np.count_nonzero(positive))


def fit_spatial_plane(
    design: np.ndarray,
    values: np.ndarray,
    covariance: np.ndarray,
    method: str,
) -> FitResult:
    if method not in METHODS:
        raise ValueError(f"unknown method {method}")
    if design.ndim != 2 or design.shape[1] != 3 or values.shape != (design.shape[0],):
        raise SpatialInputError("invalid spatial design dimensions")
    if covariance.shape != (design.shape[0], design.shape[0]):
        raise SpatialInputError("spatial covariance dimension mismatch")
    if design.shape[0] < 4:
        raise SpatialInputError("fewer than four spatial contrasts")

    if method == "none":
        beta, _, rank, _ = np.linalg.lstsq(design, values, rcond=1e-12)
        residual = values - design @ beta
        degrees_of_freedom = max(design.shape[0] - int(rank), 1)
        residual_variance = float(residual @ residual / degrees_of_freedom)
        information_inverse, information_rank = _symmetric_pseudoinverse(
            design.T @ design
        )
        beta_covariance = residual_variance * information_inverse
        return FitResult(
            method=method,
            beta=beta,
            beta_covariance=beta_covariance,
            rank=information_rank,
            training_count=design.shape[0],
            residual_variance=residual_variance,
            covariance_rank=design.shape[0],
        )

    if method == "diagonal":
        variances = np.diag(covariance)
        scale = max(1.0, float(np.max(np.abs(variances))))
        tolerance = REGULARISATION_ABSOLUTE + REGULARISATION_RELATIVE * scale
        if np.any(variances <= tolerance):
            raise SpatialInputError("nonpositive diagonal SD variance")
        covariance_inverse = np.diag(1 / variances)
        covariance_rank = design.shape[0]
    else:
        covariance_inverse, covariance_rank = _symmetric_pseudoinverse(covariance)

    information = design.T @ covariance_inverse @ design
    information_inverse, information_rank = _symmetric_pseudoinverse(information)
    beta = information_inverse @ design.T @ covariance_inverse @ values
    residual = values - design @ beta
    weighted_residual = float(residual @ covariance_inverse @ residual)
    degrees_of_freedom = max(covariance_rank - information_rank, 1)
    return FitResult(
        method=method,
        beta=beta,
        beta_covariance=information_inverse,
        rank=information_rank,
        training_count=design.shape[0],
        residual_variance=weighted_residual / degrees_of_freedom,
        covariance_rank=covariance_rank,
    )


def predict_spatial_plane(
    fit: FitResult, design: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    mean = design @ fit.beta
    variance = np.einsum("ij,jk,ik->i", design, fit.beta_covariance, design)
    variance = np.maximum(variance, 0)
    return mean, variance


def product_gap_keys(product_audit_path: Path, label: str) -> set[tuple[int, Decimal]]:
    payload = json.loads(product_audit_path.read_text(encoding="utf-8"))
    try:
        experiment = payload["experiments"][label]
        defective = experiment["validation"]["eligible_satellites_with_coverage_defects"]
    except (KeyError, TypeError) as error:
        raise SpatialInputError("unrecognised product audit structure") from error
    keys: set[tuple[int, Decimal]] = set()
    for satellite in defective:
        for gap in experiment["clk"]["satellites"][satellite]["gaps"]:
            current = datetime.fromisoformat(gap["first_epoch"])
            final = datetime.fromisoformat(gap["last_epoch"])
            while current <= final:
                seconds = Decimal(str((current - GPS_EPOCH).total_seconds()))
                week = int(seconds // Decimal(604800))
                tow = seconds - Decimal(week * 604800)
                keys.add((week, tow))
                current += timedelta(seconds=30)
    return keys


def _valid_spatial_epoch(epoch: SpatialEpoch, gaps: set[tuple[int, Decimal]]) -> bool:
    return (
        epoch.key not in gaps
        and epoch.status == "OK"
        and epoch.declared_difference_count > 0
        and epoch.dropped_geometry_count == 0
        and epoch.values.size == epoch.declared_difference_count
        and np.isfinite(epoch.design).all()
        and np.isfinite(epoch.values).all()
        and np.isfinite(epoch.covariance).all()
    )


def collect_loso_calibration(
    stec_path: Path,
    sd_path: Path,
    product_audit_path: Path,
    expected_sites: set[str],
) -> tuple[dict[str, list[CalibrationRecord]], dict[str, object]]:
    gaps = product_gap_keys(product_audit_path, "quiet")
    records: dict[str, list[CalibrationRecord]] = {method: [] for method in METHODS}
    observed_sites: set[str] = set()
    epoch_count = 0
    usable_epoch_count = 0
    dropped_geometry_epochs = 0
    fit_failures: Counter[str] = Counter()
    for epoch in iter_spatial_epochs(stec_path, sd_path):
        epoch_count += 1
        observed_sites.update(str(site) for site in epoch.sites)
        if epoch.dropped_geometry_count:
            dropped_geometry_epochs += 1
        if not _valid_spatial_epoch(epoch, gaps):
            continue
        usable_epoch_count += 1
        for site in sorted(set(str(value) for value in epoch.sites)):
            test = epoch.sites == site
            train = ~test
            if np.count_nonzero(test) == 0 or len(set(epoch.sites[train])) < 3:
                continue
            train_covariance = epoch.covariance[np.ix_(train, train)]
            target_variance = np.diag(epoch.covariance)[test]
            for method in METHODS:
                try:
                    fit = fit_spatial_plane(
                        epoch.design[train],
                        epoch.values[train],
                        train_covariance,
                        method,
                    )
                    prediction, prediction_variance = predict_spatial_plane(
                        fit, epoch.design[test]
                    )
                except SpatialInputError:
                    fit_failures[method] += 1
                    continue
                for residual, model_variance, validation_variance in zip(
                    epoch.values[test] - prediction,
                    prediction_variance,
                    target_variance,
                    strict=True,
                ):
                    records[method].append(
                        CalibrationRecord(
                            gps_week=epoch.gps_week,
                            gps_tow=epoch.gps_tow,
                            site=site,
                            residual=float(residual),
                            prediction_variance=float(model_variance),
                            target_variance=float(validation_variance),
                        )
                    )

    if observed_sites != expected_sites:
        raise SpatialInputError(
            f"quiet model sites mismatch: observed={sorted(observed_sites)}, "
            f"expected={sorted(expected_sites)}"
        )
    return records, {
        "epoch_count": epoch_count,
        "usable_epoch_count": usable_epoch_count,
        "registered_product_gap_epoch_count": len(gaps),
        "dropped_geometry_epoch_count": dropped_geometry_epochs,
        "fit_failures": dict(fit_failures),
        "observed_sites": sorted(observed_sites),
    }


def _block_equal_weights(records: list[CalibrationRecord]) -> np.ndarray:
    counts = Counter(record.block for record in records)
    weights = np.array([1 / counts[record.block] for record in records], dtype=float)
    return weights / weights.sum()


def calibration_nll(records: list[CalibrationRecord], kappa: float) -> float:
    if not records or kappa < 0 or not math.isfinite(kappa):
        return math.inf
    residual = np.array([record.residual for record in records])
    prediction = np.array([record.prediction_variance for record in records])
    target = np.array([record.target_variance for record in records])
    variance = kappa * kappa * prediction + target
    if np.any(~np.isfinite(variance)) or np.any(variance <= 0):
        return math.inf
    weights = _block_equal_weights(records)
    terms = 0.5 * (np.log(2 * math.pi * variance) + residual * residual / variance)
    return float(weights @ terms)


def optimise_kappa(records: list[CalibrationRecord]) -> tuple[float, float]:
    candidates = np.concatenate(([0.0], np.geomspace(1e-4, 1e4, 321)))
    objective = np.array([calibration_nll(records, float(value)) for value in candidates])
    best = int(np.argmin(objective))
    if best == 0 or best == len(candidates) - 1:
        return float(candidates[best]), float(objective[best])
    low = math.log(candidates[best - 1])
    high = math.log(candidates[best + 1])
    golden = (math.sqrt(5) - 1) / 2
    left = high - golden * (high - low)
    right = low + golden * (high - low)
    for _ in range(60):
        if calibration_nll(records, math.exp(left)) < calibration_nll(
            records, math.exp(right)
        ):
            high = right
            right = left
            left = high - golden * (high - low)
        else:
            low = left
            left = right
            right = low + golden * (high - low)
    kappa = math.exp((low + high) / 2)
    return kappa, calibration_nll(records, kappa)


def calibration_summary(
    records: list[CalibrationRecord], kappa: float
) -> dict[str, object]:
    residual = np.array([record.residual for record in records])
    prediction = np.array([record.prediction_variance for record in records])
    target = np.array([record.target_variance for record in records])
    variance = kappa * kappa * prediction + target
    half_width = 1.959963984540054 * np.sqrt(variance)
    coverage = float(np.mean(np.abs(residual) <= half_width))
    return {
        "record_count": len(records),
        "epoch_count": len({(record.gps_week, record.gps_tow) for record in records}),
        "site_count": len({record.site for record in records}),
        "one_hour_block_count": len({record.block for record in records}),
        "kappa": kappa,
        "block_equal_weight_nll": calibration_nll(records, kappa),
        "pointwise_95_coverage": coverage,
        "mean_interval_width_tecu": float(np.mean(2 * half_width)),
        "rmse_tecu": float(np.sqrt(np.mean(residual * residual))),
    }


def build_freeze(
    quiet_stec: Path,
    quiet_sd: Path,
    product_audit: Path,
    validation_split_audit: Path,
) -> dict[str, object]:
    split = json.loads(validation_split_audit.read_text(encoding="utf-8"))
    if not split.get("valid") or not split.get("model_split_valid"):
        raise SpatialInputError("validation split audit is not valid")
    expected_sites = set(split["model_stations"])
    records, source_summary = collect_loso_calibration(
        quiet_stec, quiet_sd, product_audit, expected_sites
    )
    calibration: dict[str, object] = {}
    for method in METHODS:
        if not records[method]:
            raise SpatialInputError(f"no LOSO calibration records for {method}")
        kappa, _ = optimise_kappa(records[method])
        calibration[method] = calibration_summary(records[method], kappa)

    return {
        "schema": FREEZE_SCHEMA,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "heldout_output_read": False,
        "scientific_scope": "FLOAT_STEC_COVARIANCE_PIPELINE_PIVOT_NOT_FIXED_STEC",
        "software": {
            "script": str(Path(__file__).resolve()),
            "methods": list(METHODS),
        },
        "geometry": {
            "earth_radius_m": EARTH_RADIUS_M,
            "shell_height_m": SHELL_HEIGHT_M,
            "mapping_function": "MSLM",
            "mslm_alpha": MSLM_ALPHA,
            "minimum_elevation_deg": math.degrees(MIN_ELEVATION_RAD),
            "origin_latitude_deg": math.degrees(ORIGIN_LAT_RAD),
            "origin_longitude_deg": math.degrees(ORIGIN_LON_RAD),
            "basis": ["intercept", "east_1000km", "north_1000km"],
            "response": "satellite_minus_reference_STEC_TE CU".replace(" ", ""),
            "design_row": "F_target*h(IPP_target)-F_reference*h(IPP_reference)",
        },
        "covariance_arms": {
            "none": "OLS; training SD covariance ignored",
            "diagonal": "GLS with diag(Q_SD)",
            "full": "GLS with complete Q_SD",
            "validation_target": "all arms add heldout FLOAT SD marginal variance",
        },
        "calibration_protocol": {
            "source": "quiet model stations only",
            "scheme": "leave-one-site-out at each epoch; continuous one-hour blocks receive equal NLL weight",
            "objective": "minimise NLL of residual variance kappa^2*V_prediction+Q_validation",
            "storm_recalibration": False,
        },
        "source_summary": source_summary,
        "calibration": calibration,
        "input_hashes": {
            "quiet_model_stec": {
                "path": str(quiet_stec.resolve()),
                "sha256": _sha256(quiet_stec),
            },
            "quiet_model_sd": {
                "path": str(quiet_sd.resolve()),
                "sha256": _sha256(quiet_sd),
            },
            "product_audit": {
                "path": str(product_audit.resolve()),
                "sha256": _sha256(product_audit),
            },
            "validation_split_audit": {
                "path": str(validation_split_audit.resolve()),
                "sha256": _sha256(validation_split_audit),
            },
        },
        "registered_thresholds": {
            "full_vs_none": {
                "coverage_error_reduction_points": 5.0,
                "or_interval_score_reduction_percent": 10.0,
                "maximum_width_increase_percent": 25.0,
            },
            "full_vs_diagonal_materiality": {
                "coverage_difference_points": 3.0,
                "nll_difference_percent": 5.0,
            },
            "calibrated_coverage": {
                "overall_percent": [93.0, 97.0],
                "worst_site_minimum_percent": 90.0,
                "maximum_width_inflation_percent": 40.0,
            },
            "quiet_to_storm_maximum_coverage_loss_points": 5.0,
            "interpretation": "diagnostic analogues only because G2 failed",
        },
    }


def _verified_freeze(path: Path) -> dict[str, object]:
    freeze = json.loads(path.read_text(encoding="utf-8"))
    if freeze.get("schema") != FREEZE_SCHEMA or freeze.get("heldout_output_read") is not False:
        raise SpatialInputError("invalid or contaminated spatial freeze")
    for record in freeze["input_hashes"].values():
        source = Path(record["path"])
        if _sha256(source) != record["sha256"]:
            raise SpatialInputError(f"frozen input hash changed: {source}")
    return freeze


def _paired_spatial_epochs(
    model_stec: Path,
    model_sd: Path,
    heldout_stec: Path,
    heldout_sd: Path,
) -> Iterator[tuple[SpatialEpoch, SpatialEpoch]]:
    model_iterator = iter_spatial_epochs(model_stec, model_sd)
    heldout_iterator = iter_spatial_epochs(heldout_stec, heldout_sd)
    for model, heldout in zip(model_iterator, heldout_iterator, strict=True):
        if model.key != heldout.key:
            raise SpatialInputError(
                f"model/heldout epoch mismatch: {model.key} versus {heldout.key}"
            )
        yield model, heldout


def predict_day(
    label: str,
    model_stec: Path,
    model_sd: Path,
    heldout_stec: Path,
    heldout_sd: Path,
    product_audit: Path,
    expected_model_sites: set[str],
    expected_heldout_sites: set[str],
    kappas: dict[str, float],
) -> tuple[list[PredictionRow], dict[str, object]]:
    gaps = product_gap_keys(product_audit, label)
    rows: list[PredictionRow] = []
    epoch_count = 0
    usable_epoch_count = 0
    model_sites: set[str] = set()
    heldout_sites: set[str] = set()
    skipped = Counter()
    for model, heldout in _paired_spatial_epochs(
        model_stec, model_sd, heldout_stec, heldout_sd
    ):
        epoch_count += 1
        model_sites.update(str(site) for site in model.sites)
        heldout_sites.update(str(site) for site in heldout.sites)
        if model.key in gaps:
            skipped["registered_product_gap"] += 1
            continue
        if not _valid_spatial_epoch(model, set()):
            skipped["invalid_model_epoch"] += 1
            continue
        if not _valid_spatial_epoch(heldout, set()):
            skipped["invalid_heldout_epoch"] += 1
            continue
        usable_epoch_count += 1
        fitted: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        for method in METHODS:
            fit = fit_spatial_plane(
                model.design, model.values, model.covariance, method
            )
            fitted[method] = predict_spatial_plane(fit, heldout.design)
        target_variance = np.diag(heldout.covariance)
        for index in range(heldout.values.size):
            method_payload: dict[str, dict[str, float]] = {}
            for method in METHODS:
                mean, prediction_variance = fitted[method]
                base = float(prediction_variance[index])
                target = float(target_variance[index])
                uncalibrated = base + target
                calibrated = kappas[method] ** 2 * base + target
                method_payload[method] = {
                    "predicted_tecu": float(mean[index]),
                    "prediction_variance_tecu2": base,
                    "target_variance_tecu2": target,
                    "uncalibrated_total_variance_tecu2": uncalibrated,
                    "calibrated_total_variance_tecu2": calibrated,
                }
            rows.append(
                PredictionRow(
                    gps_week=heldout.gps_week,
                    gps_tow=heldout.gps_tow,
                    site=str(heldout.sites[index]),
                    target_satellite=heldout.target_satellites[index],
                    reference_satellite=heldout.reference_satellites[index],
                    observed_tecu=float(heldout.values[index]),
                    methods=method_payload,
                )
            )

    if model_sites != expected_model_sites:
        raise SpatialInputError(
            f"{label} model sites mismatch: {sorted(model_sites)}"
        )
    if heldout_sites != expected_heldout_sites:
        raise SpatialInputError(
            f"{label} heldout sites mismatch: {sorted(heldout_sites)}"
        )
    return rows, {
        "epoch_count": epoch_count,
        "usable_epoch_count": usable_epoch_count,
        "prediction_row_count": len(rows),
        "registered_product_gap_epoch_count": len(gaps),
        "skipped_epochs": dict(skipped),
        "model_sites": sorted(model_sites),
        "heldout_sites": sorted(heldout_sites),
    }


def _metric_values(
    rows: list[PredictionRow], method: str, calibrated: bool
) -> dict[str, float | int]:
    if not rows:
        raise SpatialInputError("cannot score zero prediction rows")
    observed = np.array([row.observed_tecu for row in rows])
    predicted = np.array([row.methods[method]["predicted_tecu"] for row in rows])
    variance_field = (
        "calibrated_total_variance_tecu2"
        if calibrated
        else "uncalibrated_total_variance_tecu2"
    )
    variance = np.array([row.methods[method][variance_field] for row in rows])
    if np.any(~np.isfinite(variance)) or np.any(variance <= 0):
        raise SpatialInputError("nonpositive prediction variance")
    residual = observed - predicted
    standard = np.sqrt(variance)
    half_width = 1.959963984540054 * standard
    covered = np.abs(residual) <= half_width
    lower = predicted - half_width
    upper = predicted + half_width
    alpha = 0.05
    interval_score = upper - lower
    interval_score += np.where(observed < lower, 2 / alpha * (lower - observed), 0)
    interval_score += np.where(observed > upper, 2 / alpha * (observed - upper), 0)
    nll = 0.5 * (np.log(2 * math.pi * variance) + residual * residual / variance)
    return {
        "row_count": len(rows),
        "epoch_count": len({(row.gps_week, row.gps_tow) for row in rows}),
        "bias_tecu": float(np.mean(residual)),
        "rmse_tecu": float(np.sqrt(np.mean(residual * residual))),
        "nll": float(np.mean(nll)),
        "coverage_percent": float(100 * np.mean(covered)),
        "mean_width_tecu": float(np.mean(2 * half_width)),
        "mean_interval_score_tecu": float(np.mean(interval_score)),
    }


def score_rows(rows: list[PredictionRow]) -> dict[str, object]:
    output: dict[str, object] = {}
    sites = sorted({row.site for row in rows})
    for method in METHODS:
        per_site = {
            site: {
                "uncalibrated": _metric_values(
                    [row for row in rows if row.site == site], method, False
                ),
                "calibrated": _metric_values(
                    [row for row in rows if row.site == site], method, True
                ),
            }
            for site in sites
        }
        overall_uncalibrated = _metric_values(rows, method, False)
        overall_calibrated = _metric_values(rows, method, True)
        output[method] = {
            "uncalibrated": overall_uncalibrated,
            "calibrated": overall_calibrated,
            "width_inflation_percent": 100
            * (
                overall_calibrated["mean_width_tecu"]
                / overall_uncalibrated["mean_width_tecu"]
                - 1
            ),
            "worst_site_calibrated_coverage_percent": min(
                payload["calibrated"]["coverage_percent"]
                for payload in per_site.values()
            ),
            "per_site": per_site,
        }
    return output


def comparison_metrics(scores: dict[str, object]) -> dict[str, float | bool]:
    none = scores["none"]["calibrated"]
    diagonal = scores["diagonal"]["calibrated"]
    full = scores["full"]["calibrated"]
    coverage_error_reduction = (
        abs(none["coverage_percent"] - 95) - abs(full["coverage_percent"] - 95)
    )
    interval_score_reduction = 100 * (
        none["mean_interval_score_tecu"] - full["mean_interval_score_tecu"]
    ) / none["mean_interval_score_tecu"]
    width_increase = 100 * (
        full["mean_width_tecu"] / none["mean_width_tecu"] - 1
    )
    coverage_difference = abs(
        full["coverage_percent"] - diagonal["coverage_percent"]
    )
    nll_difference = 100 * (diagonal["nll"] - full["nll"]) / max(
        abs(diagonal["nll"]), 1e-12
    )
    return {
        "full_vs_none_coverage_error_reduction_points": coverage_error_reduction,
        "full_vs_none_interval_score_reduction_percent": interval_score_reduction,
        "full_vs_none_width_increase_percent": width_increase,
        "full_vs_none_threshold_passed": (
            (coverage_error_reduction >= 5 or interval_score_reduction >= 10)
            and width_increase <= 25
        ),
        "full_vs_diagonal_coverage_difference_points": coverage_difference,
        "full_vs_diagonal_nll_improvement_percent": nll_difference,
        "full_vs_diagonal_material_benefit": (
            coverage_difference >= 3 or nll_difference >= 5
        ),
    }


def _bootstrap_comparisons(rows: list[PredictionRow]) -> dict[str, object]:
    by_block: dict[str, list[PredictionRow]] = defaultdict(list)
    for row in rows:
        by_block[row.block].append(row)
    blocks = sorted(by_block)
    generator = np.random.default_rng(BOOTSTRAP_SEED)
    samples: dict[str, list[float]] = defaultdict(list)
    for _ in range(BOOTSTRAP_REPLICATES):
        selected = generator.choice(blocks, size=len(blocks), replace=True)
        sample = [row for block in selected for row in by_block[str(block)]]
        comparison = comparison_metrics(score_rows(sample))
        for key, value in comparison.items():
            if isinstance(value, bool):
                continue
            samples[key].append(float(value))
    return {
        "scheme": "resample contiguous one-hour epoch blocks with replacement",
        "replicates": BOOTSTRAP_REPLICATES,
        "seed": BOOTSTRAP_SEED,
        "block_count": len(blocks),
        "percentile_95_intervals": {
            key: {
                "lower": float(np.percentile(values, 2.5)),
                "upper": float(np.percentile(values, 97.5)),
            }
            for key, values in samples.items()
        },
    }


def build_report(
    freeze_path: Path,
    product_audit: Path,
    quiet_model_stec: Path,
    quiet_model_sd: Path,
    quiet_heldout_stec: Path,
    quiet_heldout_sd: Path,
    storm_model_stec: Path,
    storm_model_sd: Path,
    storm_heldout_stec: Path,
    storm_heldout_sd: Path,
) -> dict[str, object]:
    freeze = _verified_freeze(freeze_path)
    if _sha256(product_audit) != freeze["input_hashes"]["product_audit"]["sha256"]:
        raise SpatialInputError("product audit differs from frozen copy")
    if _sha256(quiet_model_stec) != freeze["input_hashes"]["quiet_model_stec"]["sha256"]:
        raise SpatialInputError("quiet model STEC differs from freeze")
    if _sha256(quiet_model_sd) != freeze["input_hashes"]["quiet_model_sd"]["sha256"]:
        raise SpatialInputError("quiet model SD differs from freeze")

    split = json.loads(
        Path(freeze["input_hashes"]["validation_split_audit"]["path"]).read_text(
            encoding="utf-8"
        )
    )
    model_sites = set(split["model_stations"])
    heldout_sites = set(split["heldout_stations"])
    kappas = {
        method: float(freeze["calibration"][method]["kappa"])
        for method in METHODS
    }
    paths = {
        "quiet": (
            quiet_model_stec,
            quiet_model_sd,
            quiet_heldout_stec,
            quiet_heldout_sd,
        ),
        "storm": (
            storm_model_stec,
            storm_model_sd,
            storm_heldout_stec,
            storm_heldout_sd,
        ),
    }
    report_days: dict[str, object] = {}
    day_scores: dict[str, dict[str, object]] = {}
    for label, day_paths in paths.items():
        rows, source = predict_day(
            label,
            *day_paths,
            product_audit,
            model_sites,
            heldout_sites,
            kappas,
        )
        scores = score_rows(rows)
        day_scores[label] = scores
        report_days[label] = {
            "source": source,
            "scores": scores,
            "comparisons": comparison_metrics(scores),
            "block_bootstrap": _bootstrap_comparisons(rows),
        }

    quiet_full = day_scores["quiet"]["full"]["calibrated"]
    storm_full = day_scores["storm"]["full"]["calibrated"]
    coverage_loss = quiet_full["coverage_percent"] - storm_full["coverage_percent"]
    return {
        "schema": REPORT_SCHEMA,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scientific_scope": freeze["scientific_scope"],
        "freeze": {
            "path": str(freeze_path.resolve()),
            "sha256": _sha256(freeze_path),
            "heldout_output_read_when_frozen": freeze["heldout_output_read"],
        },
        "days": report_days,
        "transfer": {
            "quiet_to_storm_full_calibrated_coverage_loss_points": coverage_loss,
            "requires_storm_conditioning": coverage_loss > 5,
        },
        "limitations": [
            "G2 failed on both AR replications; this report does not test fixed STEC.",
            "Held-out FLOAT estimates are noisy validation targets, not independent geophysical truth.",
            "The fitted east/north plane coefficients are model-dependent gradient proxies.",
            "No held-out observation, covariance, residual or score was used to choose the model or kappa.",
        ],
        "input_hashes": {
            name: {"path": str(path.resolve()), "sha256": _sha256(path)}
            for name, path in {
                "product_audit": product_audit,
                "quiet_model_stec": quiet_model_stec,
                "quiet_model_sd": quiet_model_sd,
                "quiet_heldout_stec": quiet_heldout_stec,
                "quiet_heldout_sd": quiet_heldout_sd,
                "storm_model_stec": storm_model_stec,
                "storm_model_sd": storm_model_sd,
                "storm_heldout_stec": storm_heldout_stec,
                "storm_heldout_sd": storm_heldout_sd,
            }.items()
        },
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    freeze_parser = subparsers.add_parser("freeze")
    freeze_parser.add_argument("--quiet-stec", type=Path, required=True)
    freeze_parser.add_argument("--quiet-sd", type=Path, required=True)
    freeze_parser.add_argument("--product-audit", type=Path, required=True)
    freeze_parser.add_argument("--validation-split-audit", type=Path, required=True)
    freeze_parser.add_argument("--output", type=Path, required=True)

    evaluate_parser = subparsers.add_parser("evaluate")
    evaluate_parser.add_argument("--freeze", type=Path, required=True)
    evaluate_parser.add_argument("--product-audit", type=Path, required=True)
    for label in ("quiet", "storm"):
        evaluate_parser.add_argument(f"--{label}-model-stec", type=Path, required=True)
        evaluate_parser.add_argument(f"--{label}-model-sd", type=Path, required=True)
        evaluate_parser.add_argument(f"--{label}-heldout-stec", type=Path, required=True)
        evaluate_parser.add_argument(f"--{label}-heldout-sd", type=Path, required=True)
    evaluate_parser.add_argument("--output", type=Path, required=True)

    args = parser.parse_args(argv)
    try:
        if args.command == "freeze":
            payload = build_freeze(
                args.quiet_stec,
                args.quiet_sd,
                args.product_audit,
                args.validation_split_audit,
            )
        else:
            payload = build_report(
                args.freeze,
                args.product_audit,
                args.quiet_model_stec,
                args.quiet_model_sd,
                args.quiet_heldout_stec,
                args.quiet_heldout_sd,
                args.storm_model_stec,
                args.storm_model_sd,
                args.storm_heldout_stec,
                args.storm_heldout_sd,
            )
    except (SpatialInputError, OSError, UnicodeError, json.JSONDecodeError) as error:
        print(f"spatial validation input error: {error}")
        return 2
    _atomic_write_json(args.output, payload)
    print(json.dumps(payload, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
