# Canine somatic VCF and tumor RNA

Vaxrank can process a dog somatic VCF and tumor RNA BAM when both use the same
canine assembly. Supply the dog's measured DLA class I alleles and an Ensembl
annotation for that assembly. This example pins CanFam3.1 and Ensembl 98; use
an annotation matched to your own files if they use a different assembly.

Install and index the canine annotation before the run:

```sh
pyensembl install --release 98 --species dog
mhcflurry-downloads fetch
```

MHCflurry model weights are a separate download. Confirm that the installed
model supports every typed allele and peptide length before interpreting its
output. An executable DLA prediction alone does not establish canine binding
accuracy or calibration; see the [MHCflurry DLA evidence audit](https://openvax.github.io/mhcflurry/dla_validation.html).

With a matched-normal somatic VCF, tumor RNA BAM and typed DLA allele, run:

```sh
mkdir -p dog-results
vaxrank \
  --vcf dog-somatic.vcf --bam dog-tumor-rna.bam \
  --genome CanFam3.1 --ensembl-release 98 --host-species dog \
  --mhc-predictor mhcflurry --mhc-alleles 'DLA-88*034:01' \
  --output-ascii-report dog-results/report.txt \
  --output-epitopes dog-results/epitopes.csv \
  --output-isovar-csv dog-results/rna-support.csv
```

Replace the example DLA allele with the dog's actual genotype, including each
allele to evaluate. The VCF must contain the called somatic variants; Vaxrank
does not call them from DNA reads. The BAM must contain tumor RNA alignments to
the VCF's assembly. The Ensembl annotation supplies translated proteins for
reference-peptide checks; incomplete annotation limits those checks. Inspect
RNA support, allele identity and reference scope along with the epitope scores.

For `--json-variants` containing a custom PyEnsembl `Genome`, its serialized
metadata may omit species and taxonomy. Supply an explicit ORFdata context:

```python
import json
from orfdata import ReferenceContext

context = ReferenceContext(9615, 'dog', 'CanFam3.1', 'Ensembl', '98')
with open('dog-reference-context.json', 'w') as stream:
    json.dump(context.to_dict(), stream)
```

Use `--reference-context dog-reference-context.json` with `--json-variants`,
`--bam` and `--host-species dog`. Vaxrank checks the declaration against every
assembly, annotation, species and taxonomy field present in the custom Genome
before acquiring an MHC predictor. The context declares missing metadata; it
does not infer species from `CanFam3.1` or verify the source FASTA/GTF bytes.

For a small reproducible run using a published dog tumor/normal call and
original tumor RNA, see the [WLOS_0001 TP53 example](https://github.com/openvax/vaxrank/tree/main/examples/wlos0001_tp53).
Its default random predictor tests the pipeline offline; the optional
MHCflurry run exercises a DLA model with an illustrative allele. WLOS_0001's
DLA genotype is unknown and that fixture contains only three TP53 reference
proteins, so neither run is a whole-sample canine analysis.

Canine CTA expression inputs use the separate experimental
[host-reference workflow](antigen-inputs.md#species-specific-cta-expression-experimental).
