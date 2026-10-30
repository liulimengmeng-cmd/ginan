# Independent product-segment numerical regression

`four_star_cases.inc` contains eight small regression inputs derived from the
frozen C39/C36/C27/C33 product records and an independently reconstructed Decimal
oracle. `numerical_core_blackbox.cpp` tests the production numerical/publication
helpers against that oracle; the helper is not used to construct the expected
values. The regression also rejects an arbitrary shift at a very large gamma.

From the repository root on Linux:

```sh
c++ -std=c++20 -O1 -UNDEBUG -Isrc/cpp \
  src/tests/unit/product_numerical_continuity.cpp -o /tmp/product_numerical_continuity
/tmp/product_numerical_continuity
c++ -std=c++20 -O1 -UNDEBUG -Isrc/cpp \
  src/tests/regression/product_segment_repair/numerical_core_blackbox.cpp \
  -o /tmp/product_segment_oracle
/tmp/product_segment_oracle
```

Expected totals are 262 and 116 checks, respectively. These are contract tests,
not GNSS epoch replay, scientific product validation or full writer-process
failure injection. The tiny-rename boundary test is separately registered as
`zhang_product_publication_boundary_tests` in CMake.
