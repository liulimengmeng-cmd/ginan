"""Syntax checks or real PEA-linked covariance-reader contract; no solver run."""
from pathlib import Path
import argparse
import os
import shlex
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("build", type=Path)
parser.add_argument("mode", choices=("syntax", "covariance"))
args = parser.parse_args()
root = Path(__file__).resolve().parents[2]
build = args.build.resolve()
output = build / "multi_gnss_clock_contract"
output.mkdir(exist_ok=True)
flags = {}
for line in (build / "cpp/CMakeFiles/pea.dir/flags.make").read_text().splitlines():
    if " = " in line:
        key, value = line.split(" = ", 1)
        flags[key] = shlex.split(value)
link = shlex.split((build / "cpp/CMakeFiles/pea.dir/link.txt").read_text())
prefix = [link[0], *flags["CXX_DEFINES"], *flags["CXX_INCLUDES"],
          *flags["CXX_FLAGS"], "-O0", "-std=c++20"]
if args.mode == "syntax":
    for name in ("ppp_obs.cpp", "ppppp.cpp", "ppp.cpp", "zhangPppAr.cpp"):
        logfile = output / (name + ".syntax.log")
        command = [*prefix, "-fsyntax-only", str(root / "src/cpp/pea" / name)]
        with logfile.open("w") as stream:
            result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT)
        print(name, "PASS" if result.returncode == 0 else "FAIL", logfile, flush=True)
        if result.returncode:
            print(logfile.read_text()[-12000:], flush=True)
            raise SystemExit(result.returncode)
else:
    driver = root / "src/tests/unit/multi_gnss_clock_product_covariance.cpp"
    subprocess.run([*prefix, "-c", str(driver), "-o", str(output / "driver.o")], check=True)
    subprocess.run(["objcopy", "--redefine-sym", "main=pea_original_main",
        str(build / "cpp/CMakeFiles/pea.dir/pea/main.cpp.o"),
        str(output / "original_main.o")], check=True)
    link[link.index("CMakeFiles/pea.dir/pea/main.cpp.o")] = str(output / "original_main.o")
    # Call the actual private product writer in the contract without adding a
    # public production test hook. Only a copied object has its local symbol
    # exported/renamed; the real PEA object and binary remain unchanged.
    writer_object = build / "cpp/CMakeFiles/pea.dir/pea/zhangPppAr.cpp.o"
    symbols = subprocess.check_output(["nm", str(writer_object)], text=True)
    writer_symbols = [line.split()[-1] for line in symbols.splitlines()
                      if "appendProductERK20ZhangInternalProduct" in line
                      and ".cold" not in line and ".constprop" not in line]
    if len(writer_symbols) != 1:
        raise SystemExit("cannot uniquely bind the actual private appendProduct writer")
    writer_symbol = writer_symbols[0]
    copied_writer = output / "product_writer.o"
    subprocess.run(["objcopy", "--globalize-symbol", writer_symbol,
        str(writer_object), str(copied_writer)], check=True)
    subprocess.run(["objcopy", "--redefine-sym",
        writer_symbol + "=zhang_multi_gnss_test_append_product", str(copied_writer)], check=True)
    link[link.index("CMakeFiles/pea.dir/pea/zhangPppAr.cpp.o")] = str(copied_writer)
    link.insert(1, str(output / "driver.o"))
    link[link.index("-o") + 1] = str(output / "product_covariance_contract")
    subprocess.run(link, cwd=build / "cpp", check=True)
    subprocess.run([str(output / "product_covariance_contract"), str(output / "fixtures")],
        check=True, env={**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"})
