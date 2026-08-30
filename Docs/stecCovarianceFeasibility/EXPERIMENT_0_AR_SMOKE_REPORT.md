# Experiment 0 AR smoke report

## Result

Main's standard `LAMBDA_ALT` ambiguity-resolution routine was invoked on 60 forward-filter epochs
and 60 RTS-smoothed epochs. Every epoch had eligible ambiguity states, but the algorithm resolved
zero integer combinations and submitted zero ambiguity pseudo-observations. These outputs are not
fixed STEC.

All 120 associated `IONO_STEC` covariance blocks passed the independent structural, finite-value,
diagonal-consistency and positive-semidefinite checks. The negative AR result is therefore not a
covariance-file failure.

## Software and configuration

- upstream base: `origin/main` at `7baa32a`;
- branch: `codex/main-stec-covariance-feasibility`;
- configuration: `experiment_0_ar_smoke.yaml`;
- observations: ALIC, DARW and HOB2, 2019 day 199, GPS only;
- epochs: 60 at 30 seconds, from GPS TOW 345600 through 347370;
- AR mode: `LAMBDA_ALT`, once per epoch, no fix-and-hold;
- configured bootstrap success threshold: 0.9999;
- configured solution-ratio threshold: 3.

The run used the external product files already present under `/home/rx/GINAN/inputData`, including
`IGS2R03FIN_20191990000_01D_01D_OSB.BIA`. Product presence and successful parsing do not by
themselves prove receiver-side integer compatibility or correct ambiguity fixing.

## Forward-filter diagnosis

| Diagnostic status | Epoch count |
|---|---:|
| `RATIO_BELOW_THRESHOLD` | 56 |
| `SUCCESS_RATE_BELOW_THRESHOLD` | 3 |
| `INSUFFICIENT_DECORRELATED_AMBIGUITIES` | 1 |

- eligible ambiguities per epoch: 54–60;
- selected decorrelated ambiguities: 1–48;
- bootstrap success: 0.955644723113085–0.999998260174229;
- actual solution ratio for the 56 epochs reaching the ratio test:
  1.00073181935735–1.10868795299324;
- resolved-combination epochs: 0/60;
- minimum STEC-covariance eigenvalue across the run: 0.0177066074649684 TECU².

## RTS-smoothed diagnosis

All 60 smoothed epochs reached the candidate ratio test and were rejected as
`RATIO_BELOW_THRESHOLD`.

- eligible ambiguities per epoch: 54–60;
- selected decorrelated ambiguities: 47–51;
- bootstrap success: 0.999961860227033–0.999977071416828;
- actual solution ratio: 1.00005562070627–1.00696146563823;
- resolved-combination epochs: 0/60;
- minimum STEC-covariance eigenvalue across the run: 0.0176521746782071 TECU².

The RTS bootstrap success exceeded the configured threshold, but the candidate separation was far
below the configured ratio threshold of 3. Lowering the threshold merely to force a fix would not
be scientific validation and is not an accepted next step.

## Artifact checksums

Runtime result root: `D:\tec\exp0-ar-smoke-v2-diagnostics`.

| Artifact | Bytes | SHA-256 |
|---|---:|---|
| forward `.STEC.COV` | 916,053 | `6163362A45F57CEC3A55ABE7320712FB254497451443DF66D1830A79BBAB18C7` |
| forward validation JSON | 76,388 | `01C0A6BFCEFED007BFE77CA8594FCFA8DE377434684893C2023E74049A57D893` |
| RTS `.STEC_smoothed.COV` | 921,690 | `8F804A0ABA1A7551B15A44A25E74F870AB8DD33C59DC7ED108682613440BCBF2` |
| RTS validation JSON | 76,564 | `C5B094D45AC296A71B5B1BADA14E7C566A6B3B9D181288A3F983F4057DFF9CF4` |
| AR configuration | 1,228 | `E17947CCD403B85A4223DC7BDD25ED4AC3A523DEC121B16B007692C499A464DA` |
| test PEA binary | 14,176,840 | `BD5B53D759FC215774620611534494F3E8CC2C4284C95D18D107D87749A18CC6` |

## Consequence for the research plan

The full STEC posterior covariance prerequisite is now operational on main. The fixed-STEC
prerequisite is not. A fixed-versus-float covariance comparison and regional-gradient prediction
interval experiment would be invalid with this run because both nominal branches remained float.

The next investigation must explain the weak candidate separation without relaxing validation
until fixes appear. Required checks include the ambiguity observable definition, external phase-OSB
application and sign, receiver ambiguity datum, code/phase signal mapping, arc resets and an
independent restart or held-out correctness test.
