"""Compile/run focused pure publication-gate contracts; never links or runs PEA.

The output directory must be new. Small source identities are recorded before
and after both tests, so a concurrent edit cannot be reported as tested.
"""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time


def identity(path):
    data = path.read_bytes()
    return {"path": str(path), "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--compiler", default="/usr/bin/g++")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    sources = [
        Path(__file__).resolve(),
        root / "src/cpp/common/zhangProductNumericalContinuity.hpp",
        root / "src/cpp/common/zhangProductPublicationGate.hpp",
        root / "src/cpp/common/zhangPublishedContinuity.hpp",
        root / "src/cpp/pea/zhangPppAr.cpp",
        root / "src/tests/unit/product_numerical_continuity.cpp",
        root / "src/tests/unit/user_published_continuity_contract.cpp",
    ]
    before = [identity(path) for path in sources]
    started = time.monotonic()
    report = {"schema": "PRODUCT_NUMERICAL_CONTINUITY_TEST_V1", "passed": False,
              "scope": "Pure production header contracts; no PEA build, solver, or scientific data",
              "source_bindings_before": before, "tests": []}
    try:
        # Verify the real writer calls the tested implementation and exposes
        # unknown history in TRACE, rather than silently retaining old zeroes.
        writer = (root / "src/cpp/pea/zhangPppAr.cpp").read_text()
        required = [
            '#include "common/zhangProductNumericalContinuity.hpp"',
            "zhangProductNumericalReference(input)",
            "zhangAssessProductNumericalContinuity(",
            "zhangProductPublicationWithNumericalAssessment(",
            '" numerical_assessment="', '" numerical_assessment_reason="',
            '" has_delta="', '" history_gap_s="', '" history_max_gap_s="',
            '" previous_epoch_gpst_seconds="', '" comparison_delta_m="',
            '"numerical_reference_gap_unassessed_new_segment"',
        ]
        missing = [token for token in required if token not in writer]
        if missing:
            raise RuntimeError("Writer integration contract missing: " + repr(missing))
        report["writer_integration_markers_present"] = True
        for name in ("product_numerical_continuity", "user_published_continuity_contract"):
            executable = output / name
            command = [args.compiler, "-std=c++20", "-O1", "-Wall", "-Wextra", "-Werror",
                       "-I", str(root / "src/cpp"),
                       str(root / "src/tests/unit" / (name + ".cpp")), "-o", str(executable)]
            with (output / (name + ".build.log")).open("w") as log:
                subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                               check=True, timeout=60)
            native = subprocess.run(["/usr/bin/prlimit", "--as=1073741824", "--cpu=30", "--",
                                     str(executable)], stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, text=True, timeout=35)
            (output / (name + ".runtime.log")).write_text(native.stdout)
            passed = native.returncode == 0 and (
                "PRODUCT_NUMERICAL_CONTINUITY_PASS" in native.stdout
                if name == "product_numerical_continuity" else
                "PASS user_published_continuity_contract" in native.stdout)
            report["tests"].append({"name": name, "passed": passed,
                "exit_code": native.returncode, "command": command,
                "native_output": native.stdout, "binary": identity(executable)})
            if not passed:
                raise RuntimeError("Focused native contract failed: " + name)
        after = [identity(path) for path in sources]
        report["source_bindings_after"] = after
        report["sources_unchanged"] = before == after
        if before != after:
            raise RuntimeError("Source changed during focused tests; results are not hash-bound")
        report["passed"] = True
    except Exception as error:
        report["error"] = type(error).__name__ + ": " + str(error)
        report["source_bindings_after"] = [identity(path) for path in sources]
    finally:
        report["elapsed_seconds"] = time.monotonic() - started
        (output / "result.json").write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
