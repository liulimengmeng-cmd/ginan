# R48 user model replay while R49 continues

Branch is an isolated client worktree based on 793814c. The server process and
its frozen binary/configuration are not changed. Source conversation: OSB,
分析解算速度并测试, conversation 6aa68697-6fac-83e9-be18-acaf20169431,
latest user-model answer 42ee9bb3-8212-43dc-b6db-3a194b280b38. The private
conversation export remains in the local audit directory.

The final experiment uses revision2/experiment.json and revision2/config.
The root experiment/config files preserve the initial failed attempts.
Raw observations, server products, binaries, build outputs and measurement
dump payloads are not included in Git; manifests identify the original files.

Runtime control:
- FLOAT: AR OFF, no feedback.
- ROUND shadow/feedback: ZHANG_USER_NAMED_ILS must be absent. The legacy user
  named-target path uses ROUND even when YAML specifies LAMBDA.
- ILS shadow/feedback: ZHANG_USER_NAMED_ILS=1 explicitly selects LAMBDA_ALT,
  uses full covariance, ratio 3 and success 0.999 preselection. Original final
  named perr and NIS gates remain. Mixed rows require exact named recovery.
- ZHANG_USER_MEASUREMENT_DUMP is optional numerical-audit output, unset in the
  full suite. It records final measurement matrices in column-major doubles.

run_suite.py [audit_root] validates hashes, refuses pre-existing output roots,
and runs serially with one BLAS/OMP thread at nice 15. Frozen source and binary
hashes are in frozen_user.json. To reproduce, provision referenced inputs and
use new experiment/output identifiers; never overwrite preserved experiments.

Repairs include Joseph covariance updates in matched configurations, rejected
ILS row cleanup, and a user-only accepted conditional-moment factor commit.
The latter records F=I-KA, b=Kz, Q=0 before publishing the posterior, labeled
USER_INTEGER_CONDITION. It is not an extra independent observation. Formal
covariance and committed constraints alone do not prove correct integers.

MW physical accumulation/eligibility separation, per-component references and
cross-epoch product covariance remain limitations. Read the final local report
and result tables for actual completion, WL/L1 submissions and paired errors.
