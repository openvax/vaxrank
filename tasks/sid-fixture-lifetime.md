# Preserve the lifetime of newly materialized Sid fixtures

Deployment of merged Hitlist PR #805 stopped before publication: six failures and 62 setup errors all reference support files removed from the process's temporary Sid fixture. All 22 source support files in the reused deployment checkout have November 3 timestamps; generated read files and rewritten manifests survived. The fresh feature checkout and isolated reproduction pass. The builder's copytree uses copy2 by default and transfers old file times into a newly created temporary fixture, making it eligible for age-based cleanup.

Build support files with fresh destination timestamps while preserving their bytes and modes. Add a regression using backdated recipe support files and check the resulting fixture contents and fresh file age. Do not suppress missing-file checks, disable cleanup, or rebuild damaged evidence silently. Publish the fix through a version-bumped PR (3.51.1), then rerun the full clean-main deployment. Keep pending docs and packaging PRs above the resulting release.

Verification: regression and related Sid tests, lint, full suite, green CI, merge, ./deploy.sh from clean main, then public PyPI metadata and artifact verification. Version 3.51.0 has not been published because its deploy correctly stopped on test failures.
