# Verified multi-GNSS source snapshot — 2026-10-30

This branch publishes the accepted local source snapshot from 2026-10-26. It
includes the accumulated multi-GNSS FLOAT/AR work and the final product-segment
repair, not only the eight-file repair delta. It is not a scientific product
release or a claim that a new one-hour AR run completed.

## Product-segment repair

- Keep the proposed publication chart separate from the numerical comparison
  chart. A known old-gamma-to-zero adjustment may be removed from the comparison
  only for a NEW_SEGMENT with the appropriate old counter/datum identity and a
  proven full native FLOAT mean/covariance map.
- Require the recorded adjustment to equal the actual writer expression exactly;
  a relative tolerance must not admit arbitrary offsets at large gamma.
- Keep SAME_SEGMENT, unknown-map, real-extra-jump and double-frequency rejection
  fail-closed. A valid new segment still declares its reset; it does not acquire
  an unproved integer bridge to the old segment.
- Restore private product runtime/ledger snapshots after an exception before the
  atomic directory rename. Immediately mark a successful publish with a noexcept
  boundary so a later fatal error cannot restore memory behind a visible bundle.

The 300-second old-arc cap and the tested ratio=2, joint-NIS-disabled,
independent-per-epoch-AR policy were not relaxed. The local WSL limit remained
20 GiB. A missing ReceiverMap definition include in the lifecycle test was also
corrected; the earlier failed test receipt was preserved locally.

## Verified scope

Two affected production translation units were rebuilt and the complete 158
production objects linked. Ten regression scopes passed, including:

- Numerical continuity: 262 checks.
- Independent four-star/eight-signal CSV-derived oracle: 116 checks.
- Production-object retained-FLOAT lifecycle and complete joint handoff contracts.
- Publication-boundary helper with real tiny directory-rename fault cases:
  five cases and 21 checks.

The independent oracle is under `src/tests/regression/product_segment_repair/`.
It does not use the production helper as its reference implementation.

The corresponding local executable SHA256 was:

`6705b2e42f3840d91d8489c539984ac53c0f9d8607ebdfcdfa3775a0eda009d3`

The executable and its external libraries are not stored in this commit. Git
retains historical CRLF or normalizes text to LF; the manifest records both
accepted raw-source and LF-equivalent SHA256 values. Every selected Git blob is
checked against those two exact alternatives; no other content difference is
accepted. Other files are inherited from the branch base.

## Explicitly not established

No real event-epoch replay, new full-hour AR run, full writer-process fault
injection, nine-section checkpoint restoration or external scientific-accuracy
acceptance was performed for this repair. The earlier local run was deliberately
stopped after 39 new epochs at 2024-07-17 GPST 01:20:00. No server executable was
replaced. Old-segment continuity is not guaranteed when evidence is unavailable.

Raw observations, bulk products, checkpoints, build objects, executables and
generated experiment outputs are not newly included in this publication.
`source_manifest.json` and `validation_summary.json` describe this snapshot and
the actual validation boundary.
