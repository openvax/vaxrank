# Original Sid RNA fusion inputs

Copied unchanged from Isovar's pinned fixtures:

- `ATP5MG--KMT2A.input.json.gz`: coding-corpus, observed annotated start,
  18-aa coding hypothesis with early stop; noncoding transcript alternatives
  remain. Two original Personalis RNA templates, not all caller support.
- `TPST1--CRCP-T1.input.json.gz`: corpus, direct ONT RNA join without a
  justified annotated frame. No protein may be invented.

The gzip files include original SAM records, annotation models, explicit
RNA-to-genome maps, exact input URLs and caller/assembly provenance. Test
checksums pin both files; tests reconstruct with Isovar rather than supplying
manually blessed result dictionaries. See Isovar's `tests/data/fusions` and
https://osteosarc.com/fusions/ for acquisition details. These examples do not
establish tumor specificity, protein expression or peptide immunogenicity.
