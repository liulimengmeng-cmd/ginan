# R48 independent user PPP-AR prefix experiment

Five stations (MARS, DYNG, NICO, BREW, JPLM) are held out from the actual 180-station R48 network input. Each station estimates its position freely using the same PRODUCT_FIXED products and full product covariance in paired FLOAT and AR cases. The existing CANONICAL_USER_IF_WL_L1 path and fixed acceptance thresholds are used. No server algorithm or running binary is modified.

The frozen prefix is 2024-07-17 00:03:00 through 00:24:30 inclusive (44 epochs). prepare.py requires the corresponding complete server checkpoint, validates every covariance triangle, checks station independence and ancillary files, and refuses to overwrite the work directory. Raw input files are referenced and hashed, never copied. The 5 upstream DGEMV warnings remain unresolved; results are diagnostic and cannot establish integer truth or final product acceptance.

Run prepare.py with OPENBLAS_NUM_THREADS=1. Work and products are saved in /mnt/c/Users/rx/Documents/GINAN/r48_user_work_20260913_prefix. Run individual users serially, with one BLAS/OMP thread, keeping server resources available. Changes to experiment scripts/configuration are committed and pushed to fork before execution.

Preflight correction: PRODUCT_FIXED does not exist before 00:03:00. The initial rejected preparation is retained in r48_user_work_20260913. Both controls begin at the first actual PRODUCT_FIXED epoch; no product relabelling or future product backfilling is used.
