# Vaccine report clarity (3.47.12)

Full reports now lead with candidate evidence, preserving the compact mode's
annotation-free behavior. Linked product and run appendices retain every emitted
sequence, source mapping and setting. Predictor/version attribution and a score
guide distinguish NetMHCpan Score_EL, MHCflurry presentation/processing, the
configured Vaxrank ranking score and Pepsickle cleavage annotations.

Paired sequences, exact target intervals, explicit missing/unequal comparators,
and readable four-column evidence tables replace the compressed ten-column grid.
Product headings use antigen names and keep manifest IDs; redundant gene lines
and single-window counts are omitted. Product target highlighting requires
unambiguous placement and matching masks among all linked candidates.

The report also fixes WT affinity lookup across predictor versions (#724).
Sequence selection, prediction, ranking and emitted construct bytes do not change.
Grouping equivalent source reconstructions remains separate work in #719.

Scientific definitions: [NetMHCpan output](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1a/3-Output_format.php)
and [MHCflurry outputs](https://openvax.github.io/mhcflurry/intro.html).
