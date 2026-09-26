# Shared fixture mechanics

Acquisition, original-record selection, transport, and integrity mechanics
come from Osteosarc (migrated at 0.2.3; the bundle is now generated with 0.9.0). Keep the reviewed recipe format, source/correction/reference
pins, data bytes and scientific expectations unchanged. Thin compatibility entry
points retain the existing commands; no sibling checkout is required.

Validate historical membership/multiplicity and source/header provenance, then
run the reconstruction/prediction/ranking regressions offline. This consumer
uses published Osteosarc 0.11.x and Topiary releases. New recipes use the
versioned panel/bundle APIs; historical adapters state their SAM-text fidelity.

Tracked by [Osteosarc #15](https://github.com/iskandr/osteosarc/issues/15).
The historical regeneration command is:

```sh
python examples/osteosarc_test_data/build.py --cache CACHE \
  --offline --output NEW_COHORTS.zip
```

Use a new output destination and a populated pinned cache for offline
regeneration. The existing checked-in scientific expectations remain the oracle.
For a common versioned recipe, the builder also accepts `--panel-recipe
recipe.json --panel-source SOURCE_ID=original.bam --output NEW_DIRECTORY
--offline`. Repeat `--panel-source` for each local original input (Vaxrank also
requires `--cache CACHE`). Panel output contains indexed BAMs, checksums, retained
record multiplicities, source/header identities and selection reasons. See the
[shared workflow](https://github.com/iskandr/osteosarc/blob/v0.9.0/docs/fixture-migration.md).
