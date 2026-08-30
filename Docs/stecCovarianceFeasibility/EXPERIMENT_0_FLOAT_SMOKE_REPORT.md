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
- implementation commit and PEA runtime version: `4e6c718` / `untagged-4e6c718`;
- configuration: `experiment_0_float_smoke.yaml`;
- observations: ALIC, DARW and HOB2, 2019 day 199, 30-second RINEX;
- processed systems: GPS only;
- ambiguity resolution: off;
- epochs: GPS week 2062, TOW 345600 through 345750.

The runtime artifact was rebuilt and rerun after the implementation commit. Windows Git reported
no tracked diff before the build. The CMake version step used Windows Git because the WSL Git
client cannot resolve the Windows absolute Git metadata path stored by this cross-platform
worktree; that path-resolution issue does not indicate a source diff.

## Commands

The configuration was checked before processing:

```bash
LD_LIBRARY_PATH=/home/rx/.local/boost-1.82/lib \
./bin/pea \
  -y Docs/stecCovarianceFeasibility/experiment_0_float_smoke.yaml \
  -a EXP0_DATA_ROOT:/home/rx/GINAN/inputData \
  -a EXP0_OUTPUT_ROOT:/mnt/d/tec/exp0-float-clean-4e6c718 \
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

Runtime result root: `D:\tec\exp0-float-clean-4e6c718`.

| Artifact | Bytes | SHA-256 |
|---|---:|---|
| forward `.STEC.COV` | 94,305 | `4F9268C8BFEC90A20DF1A899A6991E2F034F74093C8CC87F981E114E4444B831` |
| forward validation JSON | 7,252 | `37B2DF1BD8F11D506C3DA4EA15575B7DCD31BE5D27E71EF2B484C2CB1E4678C7` |
| RTS `.STEC_smoothed.COV` | 94,539 | `8F84904ADBD3BF2C88FF016231E310217CB66C46D3C17F07304CB332CA27EFF6` |
| RTS validation JSON | 7,309 | `6CE0A97A831E164E29610DA897BCC3DA746C5C8174096ED99C95F7355F83B5A7` |
| smoke configuration | 1,158 | `5348F874F03CE91BD0C915BD046CB08128E3B06835FFB434690BB51797EE633E` |
| test PEA binary | 14,176,840 | `10CEA79B4793A047FA756CE565E3F288DEBD9B0698DAB40F8B75AC7009F35A82` |

The two covariance files, two validation files, configuration and exact PEA binary above identify
the clean-implementation runtime evidence.

## Open gates

- The companion AR smoke test produced no resolved combinations; see
  `EXPERIMENT_0_AR_SMOKE_REPORT.md`.
- No fixed/float covariance comparison, wrong-fix control, held-out receiver/baseline test,
  disturbed-day test, gradient reconstruction, or prediction-interval calibration has been run.
- The present result therefore passes only the software, file-completeness and numerical-covariance
  prerequisites for this six-epoch float smoke case.
