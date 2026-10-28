# Cropped mutation intervals (#720)

Native replay can rerank a long source fragment into windows before, across,
and after the mutation. Windows before the mutation currently receive a start
coordinate beyond their length, while the end is clipped. Constructing their
typed antigen then rejects the reversed interval.

Clamp both ends to the window's closed coordinate range. An interval outside
the window becomes an empty interval at its nearest edge, which cannot overlap
any contained peptide. Preserve the existing strict-spanning interpretation of
zero-width deletion junctions, source provenance, scores, and window eligibility.

Verify every contained query interval for substitution, insertion-length, and
zero-width events against the original fragment; exercise the composed ranking
path with epitope evidence in windows on both sides. Then replay the real mixed
mutation/SV report, run lint and the full test suite, and release independently
after #714.
