# R51 product generation without expensive diagnostic decompositions

The current ratio-only production run is intended to supply products for a later
PPP-AR experiment. Diagnostic-only decompositions must not delay production.

`zhangExpensiveDiagnosticsEnabled()` defaults to false when
`ZHANG_R51_RATIO_ONLY=1`. For an explicit offline diagnostic run, set
`ZHANG_R51_EXPENSIVE_DIAGNOSTICS=1`; setting it to `0` disables diagnostics in
other policies too. This switch is not a numerical acceptance bypass.

Skipped before allocating or decomposing diagnostic matrices:

- AR summary ADOP and repeated whole/signal/WL posterior covariance rank spectra.
- LAMBDA candidate covariance, marginal diagnostic statistics, whitened spectrum
  and dominant-mode loadings. Tested rows, integers and count remain available.
- Non-enforcing candidate NIS assessment with no acceptance consumer.
- Failed safe-prefix conditional-innovation covariance reconstruction.
- Integer-basis condition-number SVD (before and after reduction), plus the
  signal-local LU and exact integer rank decompositions used only in TRACE.
- Component WL covariance spectrum and extra diagnostic LAMBDA reduction.
- Product gain spectra, IAR gain audit, frozen graph representation/pair audits,
  canonical numerical cross-section exports and full search-matrix snapshots.

The mandatory posterior rank used by `integerDatumComplete` is retained.
`R48_COMPONENT_DIAGNOSIS` is a misleading historical timer name: it includes
bridge search and cannot be disabled. Exact physical entailment, integer
compatibility/HNF, PSD/null-space numerical guards, actual LAMBDA reduction and
search, posterior conditioning, tree/arc transport, writer authorization and
product covariance are unchanged. NIS values produced by a required validity
assessment remain available, but never veto ratio-only acceptance by size.

The compact policy record declares disabled fields; skipped diagnostics must not
be interpreted as zero covariance, zero NIS, or a successful validation.

Verification:

- `zhang_full_rank_tests`, including the failed-prefix/rescue on/off comparison.
- `python3 src/tests/run_r51_algebra_transport.py build` links actual production
  objects and compares diagnostics on/off for LAMBDA and LAMBDA_ALT: integer
  rows, fixed integers, counts, ratio and search calls must be identical.
- Four-epoch cross-day production replay plus checkpoint restore, comparing
  FLOAT fingerprints, product rows and product covariance against speed1.
- Full-network timing must be measured separately; tiny fixtures do not prove
  real-network speedup or correctness of the fixed integers.
