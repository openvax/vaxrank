# Standard Hitlist dependency and optional local evidence

Install Hitlist with ordinary Vaxrank; remove the `hitlist` packaging extra and all extra-based installation guidance. Declare PyArrow directly because Vaxrank reads Hitlist Parquet artifacts itself. Preserve the existing supported Hitlist range unless a public runtime API requires a newer minimum.

Data availability remains a runtime concern. No data inspection should create caches, build indexes, download datasets, or turn unavailable evidence into a negative observation. Explicit evidence bundles retain validation, expression/reference consistency checks, capture, and offline replay.

Automatically consume current local MS observations for CTA expression predictions, using the public read-only Hitlist snapshot and query APIs. Explicit bundles take precedence. Missing, stale or unverifiable data leave support unknown; a completed empty query records zero observations. Retain query scope, rows, build metadata and availability status for portable native replay and pooled inputs. Local observations do not establish CTA specificity or tissue safety, and filtering/ranking remain explicit DSL choices.

Use a new feature PR based on released main, with a minor version bump for the new installation/runtime contract. Keep pending documentation and PEP 621 PRs consistent and renumber their unreleased versions above this release rather than restoring extras later.

Verification: required-dependency and real bundle tests, absent-data/no-network behavior, installed-data behavior if automatic discovery is selected, portable replay, lint, full suite, strict docs, wheel metadata, GitHub CI, merge, immediate clean-main deploy, public PyPI version and artifact-hash verification.
