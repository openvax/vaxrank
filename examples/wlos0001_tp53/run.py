"""Replay one published canine somatic call and original tumor RNA in Vaxrank."""

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory
from unittest.mock import patch

from isovar import ProteinSequenceCreator, run_isovar
from mhctools import RandomBindingPredictor
from pyensembl import Genome
from varcode import load_vcf

from vaxrank.core_logic import run_vaxrank
from vaxrank.epitope_config import EpitopeConfig
from vaxrank.input_scope import InputProvenance, InputScope, validate_host_taxon
from vaxrank.patient_info import PatientInfo
from vaxrank.report import TemplateDataCreator, make_ascii_report


DATA = Path(__file__).resolve().parent / "data"
ALLELE = "DLA-88*001:01"  # Workflow placeholder: this dog's genotype is unknown.


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _reference(directory):
    fasta = directory / "canfam3-chr5-tp53-masked.fa"
    with gzip.open(DATA / (fasta.name + ".gz"), "rb") as source, fasta.open("wb") as dest:
        shutil.copyfileobj(source, dest)
    genome = Genome(
        reference_name="CanFam3.1",
        annotation_name="Ensembl",
        annotation_version=98,
        gtf_path_or_url=str(DATA / "tp53-ensembl-98.gtf"),
        transcript_fasta_paths_or_urls=[str(DATA / "tp53-transcripts-ensembl-98.fa")],
        protein_fasta_paths_or_urls=[str(DATA / "tp53-proteins-ensembl-98.fa")],
        genome_fasta_path_or_url=str(fasta),
        cache_directory_path=str(directory / "pyensembl-cache"),
    )
    genome.index()
    return genome


def _source(path, source_format, sample_id):
    return InputProvenance(
        input_id=path.name,
        source_format=source_format,
        path=str(path),
        content_sha256=_sha256(path),
        scope=InputScope(
            patient_id="WLOS_0001",
            sample_id=sample_id,
            reference_assembly="CanFam3.1",
            annotation="Ensembl 98",
            host_taxon_id=9615,
        ),
    )


def run_example(output_dir, scratch_dir=None, *, allele=ALLELE, predictor="random"):
    """Write a source receipt and ASCII report; return the receipt for checks."""
    if predictor not in {"random", "mhcflurry"}:
        raise ValueError(f"Unknown predictor: {predictor}")
    validate_host_taxon(InputScope(host_taxon_id=9615, mhc_alleles=(allele,)),
                        None, "Canine example", (allele,))
    if predictor == "mhcflurry":
        from mhctools import MHCflurry
        mhc_predictor = MHCflurry([allele])
    else:
        mhc_predictor = RandomBindingPredictor([allele])
    epitope_config = (EpitopeConfig() if predictor == "mhcflurry" else
                      EpitopeConfig(score_expr="1.0", min_epitope_score=0))
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    vcf = DATA / "wlos0001-tp53.vcf"
    bam = DATA / "wlos0001-tp53.bam"
    with TemporaryDirectory(dir=scratch_dir) as directory:
        directory = Path(directory)
        genome = _reference(directory)
        variants = load_vcf(str(vcf), genome=genome)
        isovar_results = run_isovar(
            variants,
            str(bam),
            protein_sequence_creator=ProteinSequenceCreator(max_protein_sequences_per_variant=0),
        )
        # The only local reference proteins are the three pinned TP53 transcripts.
        with patch.dict(os.environ, {"VAXRANK_REF_PEPTIDES_DIR": str(directory / "ref-peptides")}):
            result = run_vaxrank(
                isovar_results,
                mhc_predictor,
                epitope_config=epitope_config,
            )

        counts = result.variant_counts()
        patient_counts = dict(counts)
        patient_counts["num_somatic_variants"] = patient_counts.pop("num_total_variants")
        patient = PatientInfo(
            patient_id="WLOS_0001",
            vcf_paths=[str(vcf)],
            bam_path=str(bam),
            mhc_alleles=[allele],
            inputs=[("Published matched-normal somatic call", str(vcf)),
                    ("Original tumor RNA alignments", str(bam))],
            input_provenance=[
                _source(vcf, "somatic_vcf", "WLOS_0001_T1"),
                _source(bam, "tumor_rna_bam", "WLOS_0001_T3"),
            ],
            **patient_counts,
        )
        template = TemplateDataCreator(
            ranked_variants_with_vaccine_peptides=result.ranked_vaccine_peptides,
            patient_info=patient,
            final_review="",
            reviewers="",
            args_for_report={
                "vaccine_type": ["peptide"],
                "mhc_predictor": predictor,
                "manufacturability": False,
                "wt_epitopes": False,
                "reference_scope": "TP53-only Ensembl 98 subset",
            },
            input_json_file=None,
        ).compute_template_data()
        make_ascii_report(template, str(output_dir / "report.txt"))

        supported = isovar_results[0]
        receipt = {
            "sample": "WLOS_0001",
            "reference": {"species": "Canis lupus familiaris", "assembly": "CanFam3.1",
                          "annotation": "Ensembl 98", "protein_scope": "three TP53 transcripts only"},
            "variant": {"contig": "5", "start": 32563056, "ref": "G", "alt": "A",
                        "dna_source": "published matched-normal tumor WGS call; DNA reads not included"},
            "rna_source": {"run": "SRR8691780", "sample": "WLOS_0001_T3",
                           "reads": {"ref": supported.num_ref_reads, "alt": supported.num_alt_reads,
                                     "other": supported.num_other_reads},
                           "fragments": {"ref": supported.num_ref_fragments,
                                         "alt": supported.num_alt_fragments,
                                         "other": supported.num_other_fragments}},
            "vaxrank_counts": counts,
            "ranked_windows": [{
                "amino_acids": peptide.amino_acids,
                "mutation_interval": [peptide.mutant_protein_fragment.mutant_amino_acid_start_offset,
                                      peptide.mutant_protein_fragment.mutant_amino_acid_end_offset],
            }
                               for _, peptides in result.ranked_vaccine_peptides
                               for peptide in peptides],
            "prediction": {"method": predictor, "allele": allele,
                           "patient_genotype": "unknown", "binding_claim": False},
            "limits": ["one TP53 locus, not a whole-sample benchmark",
                       "no normal or tumor DNA reads in this fixture",
                       "no full canine self proteome or validated dog MHC binding model"],
            "files_sha256": {path.name: _sha256(path) for path in sorted(DATA.iterdir()) if path.is_file()},
            "sources": {
                "study": "https://doi.org/10.1038/s42003-019-0487-2",
                "rna": "https://www.ncbi.nlm.nih.gov/sra/SRR8691780",
                "fixture": "https://github.com/openvax/isovar/tree/main/tests/data/wlos0001_tp53",
            },
        }
        (output_dir / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
        return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="wlos0001-tp53-output")
    parser.add_argument("--allele", default=ALLELE,
                        help="Illustrative DLA allele; WLOS_0001's genotype is unknown")
    parser.add_argument("--predictor", choices=("random", "mhcflurry"), default="random")
    args = parser.parse_args()
    receipt = run_example(args.output_dir, allele=args.allele, predictor=args.predictor)
    print("RNA alt reads: %d; ranked windows: %d" % (
        receipt["rna_source"]["reads"]["alt"], len(receipt["ranked_windows"])))
    print("Wrote %s/receipt.json and report.txt" % args.output_dir)


if __name__ == "__main__":
    main()
