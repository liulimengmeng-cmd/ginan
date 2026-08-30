# Experiment 0 float smoke report

## Result

The optional full `IONO_STEC` covariance sidecar was exercised through the complete PEA executable
on six 30-second GPS epochs. Both the forward-filter and RTS-smoothed sidecars contained six
complete covariance blocks. All twelve blocks passed the independent structural, finite-value,
diagonal-consistency and positive-semidefinite checks.

This is a file-plumbing and numerical-integrity result only. Ambiguity resolution was explicitly
disabled, so none of these blocks is evidence for fixed STEC or for improved regional-gradient
prediction intervals.

## Software and data

- upstream base: `origin/main` at `7baa32a`;
- branch: `codex/main-stec-covariance-feasibility`;
- PEA runtime version: `untagged-7baa32a-dirty` while the implementation was under test;
- configuration: `experiment_0_float_smoke.yaml`;
- observations: ALIC, DARW and HOB2, 2019 day 199, 30-second RINEX;
- processed systems: GPS only;
- ambiguity resolution: off;
- epochs: GPS week 2062, TOW 345600 through 345750.

The runtime version is intentionally recorded as dirty because this report precedes the branch
commit. A clean-commit rerun is required before treating the output as a frozen experiment
artifact.

## Commands

The configuration was checked before processing:

```bash
LD_LIBRARY_PATH=/home/rx/.local/boost-1.82/lib \
./bin/pea \
  -y Docs/stecCovarianceFeasibility/experiment_0_float_smoke.yaml \
  -a EXP0_DATA_ROOT:/home/rx/GINAN/inputData \
  -a EXP0_OUTPUT_ROOT:/mnt/d/tec/exp0-float-smoke-r2 \
  --dry-run
```

The actual smoke run used the same arguments without `--dry-run`. Each sidecar was then checked
with:

```powershell
D:\Anaconda\python.exe scripts\validate_stec_covariance.py `
  <SIDE_CAR_FILE> --json-output <SIDE_CAR_FILE>.validation.json
```

## Validation evidence

| Sidecar | Epochs | State counts | Stage label | Minimum eigenvalue | Maximum diagonal error |
|---|---:|---:|---|---:|---:|
| forward | 6/6 valid | 28, 29, 30, 30, 30, 30 | `FILTER_POSTERIOR_NO_EPOCH_AR` | 0.200910840145329 TECU² | 0 TECU² |
| RTS smoothed | 6/6 valid | 28, 29, 30, 30, 30, 30 | `RTS_SMOOTHED_POSTERIOR_NO_EPOCH_AR` | 0.200910840145329 TECU² | 0 TECU² |

Every writer status was `OK`; the largest captured matrix had 30 states and 465 upper-triangle
entries. Writer-reported maximum asymmetry was zero for every block.

## Artifact checksums

Runtime result root: `D:\tec\exp0-float-smoke-v2-diagnostics`.

| Artifact | Bytes | SHA-256 |
|---|---:|---|
| forward `.STEC.COV` | 94,305 | `4F9268C8BFEC90A20DF1A899A6991E2F034F74093C8CC87F981E114E4444B831` |
| forward validation JSON | 7,259 | `87AED71CE9AA2BE74811B21E567ADAD5EE420D2031E0B8CD7B14A5C531C1D893` |
| RTS `.STEC_smoothed.COV` | 94,539 | `8F84904ADBD3BF2C88FF016231E310217CB66C46D3C17F07304CB332CA27EFF6` |
| RTS validation JSON | 7,316 | `BA4A116789F5D15203A645D2E6A11C80AA0831981EA13129198913200495BE41` |
| smoke configuration | 1,158 | `5348F874F03CE91BD0C915BD046CB08128E3B06835FFB434690BB51797EE633E` |
| test PEA binary | 14,176,840 | `BD5B53D759FC215774620611534494F3E8CC2C4284C95D18D107D87749A18CC6` |

The PEA binary and configuration hashes will change when the implementation is committed or the
report is revised; the two covariance and two validation hashes identify the reported runtime
evidence.

## Open gates

- The companion AR smoke test produced no resolved combinations; see
  `EXPERIMENT_0_AR_SMOKE_REPORT.md`.
- No fixed/float covariance comparison, wrong-fix control, held-out receiver/baseline test,
  disturbed-day test, gradient reconstruction, or prediction-interval calibration has been run.
- The present result therefore passes only the software, file-completeness and numerical-covariance
  prerequisites for this six-epoch float smoke case.
