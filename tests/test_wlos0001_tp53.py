"""Offline canine tumor DNA/RNA replay through Isovar and Vaxrank."""

import json

from orfdata import ReferenceContext
import pytest
from varcode import load_vcf

from examples.wlos0001_tp53.run import DATA, _reference, run_example
from vaxrank.cli.entry_point import run_cli


def test_wlos0001_tp53_rna_reaches_ranked_vaxrank_report(tmp_path):
    output = tmp_path / "result"
    receipt = run_example(output, scratch_dir=tmp_path)

    assert receipt["variant"] == {
        "contig": "5", "start": 32563056, "ref": "G", "alt": "A",
        "dna_source": "published matched-normal tumor WGS call; DNA reads not included",
    }
    assert receipt["rna_source"]["reads"] == {"ref": 107, "alt": 297, "other": 1}
    assert receipt["rna_source"]["fragments"] == {"ref": 102, "alt": 289, "other": 1}
    assert receipt["vaxrank_counts"] == {
        "num_total_variants": 1,
        "num_coding_effect_variants": 1,
        "num_variants_with_rna_support": 1,
        "num_variants_with_vaccine_peptides": 1,
    }
    assert any(window["amino_acids"][slice(*window["mutation_interval"])] == "C"
               for window in receipt["ranked_windows"])
    assert receipt["prediction"]["patient_genotype"] == "unknown"
    assert receipt["prediction"]["binding_claim"] is False
    assert receipt["files_sha256"]["wlos0001-tp53.bam"] == (
        "f0943702c838bc8dcc2db67410f0dd34d31ab65654b63a87c1b6f7e303ac3e00")

    report = (output / "report.txt").read_text()
    assert "Patient ID: WLOS_0001" in report
    assert "sample: WLOS_0001_T1; reference: CanFam3.1; annotation: Ensembl 98" in report
    assert "sample: WLOS_0001_T3; reference: CanFam3.1; annotation: Ensembl 98" in report
    assert "Declared patient genotype: Unknown" in report
    assert "RNA reads supporting variant allele: 297" in report
    assert "TP53" in report
    assert "random" in report
    assert "Processing predictor: pepsickle" not in report
    assert (output / "receipt.json").is_file()


def test_wlos0001_rejects_human_mhc_allele_before_writing(tmp_path):
    output = tmp_path / "result"
    with pytest.raises(ValueError, match="allele .* conflicts"):
        run_example(output, allele="HLA-A*02:01", predictor="mhcflurry")
    assert not output.exists()


def test_wlos0001_serialized_custom_genome_accepts_declared_context(tmp_path, monkeypatch):
    """A direct VCF/RNA run checks the original custom Genome before prediction."""
    genome = _reference(tmp_path)
    variants = load_vcf(str(DATA / 'wlos0001-tp53.vcf'), genome=genome)
    source = tmp_path / 'variants.json'
    variants.write_json_file(source)
    context = tmp_path / 'context.json'
    context.write_text(json.dumps(ReferenceContext(
        9615, 'dog', 'CanFam3.1', 'Ensembl', '98').to_dict()))
    monkeypatch.setenv('VAXRANK_REF_PEPTIDES_DIR', str(tmp_path / 'reference-cache'))
    report = tmp_path / 'report.txt'
    arguments = ['--json-variants', str(source), '--bam', str(DATA / 'wlos0001-tp53.bam'),
                 '--reference-context', str(context), '--host-species', 'dog',
                 '--mhc-predictor', 'random', '--mhc-alleles', 'DLA-88*034:01',
                 '--output-ascii-report', str(report)]
    wrong = json.loads(context.read_text())
    wrong['annotation_version'] = '99'
    context.write_text(json.dumps(wrong))
    with pytest.raises(ValueError, match='annotation_version'):
        run_cli(arguments)
    assert not report.exists()
    wrong['annotation_version'] = '98'
    context.write_text(json.dumps(wrong))
    run_cli(arguments)
    assert 'VCF (somatic variants):' in report.read_text()
