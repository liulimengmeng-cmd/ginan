"""Read-only historical product comparison; no datum shifts or fitted offsets.
Usage: python3 compare_r51_product_epoch.py OLD_BUNDLE NEW_BUNDLE OUTPUT.json
This is a version comparison, not a strict paired runtime benchmark.
"""
from pathlib import Path
import csv
import json
import math
import sys

left, right, output = map(Path, sys.argv[1:4])


def load(directory, name, keys):
    with (directory / name).open() as stream:
        rows = list(csv.DictReader(stream))
    mapped = {tuple(row[k] for k in keys): row for row in rows}
    assert len(mapped) == len(rows), "Duplicate product identity"
    return mapped


keys = ("gpst_seconds", "solution", "satellite", "observable")
old, new = (load(p, "products.csv", keys) for p in (left, right))
common = old.keys() & new.keys()
numeric = ("clock_m", "phase_m", "correction_m", "clock_sigma_m", "phase_sigma_m",
           "correction_sigma_m", "clock_phase_covariance_m2")
maxima = {}
for field in numeric:
    values = []
    for key in common:
        a, b = float(old[key][field]), float(new[key][field])
        if math.isfinite(a) and math.isfinite(b):
            values.append(abs(a - b))
    maxima[field] = max(values, default=None)
flags = ("ar_valid", "pppar_usable", "product_state", "discontinuity_counter",
         "public_phase_segment_id", "integer_alignment_generation", "real_gauge_generation")
flag_changes = {field: sum(old[k][field] != new[k][field] for k in common) for field in flags}
covkeys = ("gpst_seconds", "solution", "row_satellite", "row_parameter", "row_observable",
           "column_satellite", "column_parameter", "column_observable")
oldcov, newcov = (load(p, "covariance.csv", covkeys) for p in (left, right))
covcommon = oldcov.keys() & newcov.keys()
covdifferences = []
for key in covcommon:
    a, b = float(oldcov[key]["covariance_m2"]), float(newcov[key]["covariance_m2"])
    if math.isfinite(a) and math.isfinite(b):
        covdifferences.append(abs(a - b))

result = dict(scope="Historical version comparison; no offsets fitted; not a runtime pairing",
              old=str(left), new=str(right), old_product_rows=len(old), new_product_rows=len(new),
              identical_product_keys=old.keys() == new.keys(), common_product_rows=len(common),
              maximum_absolute_differences=maxima, publication_flag_changes=flag_changes,
              old_covariance_rows=len(oldcov), new_covariance_rows=len(newcov),
              identical_covariance_keys=oldcov.keys() == newcov.keys(),
              compared_finite_covariance_entries=len(covdifferences),
              maximum_covariance_difference_m2=max(covdifferences, default=None))
output.write_text(json.dumps(result, indent=2))
print(json.dumps(result, indent=2))
