"""Tests for varcode>=12.0.1 CLI compatibility and opt-in --combine-cis-snvs (#898)."""

from argparse import Namespace
import pytest
from varcode import load_vcf

from vaxrank.cli.arg_parser import make_vaxrank_arg_parser
from vaxrank.cli.entry_point import maybe_combine_cis_snvs


def _write_phased_mnv_vcf(path, *, multi_sample=False):
    # BRAF V600 (chr7 on GRCh37:75, minus strand; 140453136 A>T + 140453135 C>T -> AC>TT = V600K)
    header = (
        "##fileformat=VCFv4.2\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n'
        '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n'
    )
    if multi_sample:
        cols = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n"
        rows = (
            "7\t140453135\t.\tC\tT\t.\tPASS\t.\tGT:PS:AF\t0|0:100:0.0\t0|1:100:0.4\n"
            "7\t140453136\t.\tA\tT\t.\tPASS\t.\tGT:PS:AF\t0|0:100:0.0\t0|1:100:0.4\n"
        )
    else:
        cols = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\n"
        rows = (
            "7\t140453135\t.\tC\tT\t.\tPASS\t.\tGT:PS:AF\t0|1:100:0.4\n"
            "7\t140453136\t.\tA\tT\t.\tPASS\t.\tGT:PS:AF\t0|1:100:0.4\n"
        )
    path.write_text(header + cols + rows, encoding="utf-8")


def test_vaxrank_arg_parser_version_and_combine_cis_snvs_flag(capsys):
    parser = make_vaxrank_arg_parser()
    with pytest.raises(SystemExit) as exc:
        parser.parse_args(["--version"])
    assert exc.value.code == 0
    assert capsys.readouterr().out.strip().startswith("Vaxrank ")

    base_args = ["--bam", "tumor.bam", "--mhc-predictor", "mhcflurry"]
    args = parser.parse_args(base_args)
    assert args.combine_cis_snvs is False
    args_on = parser.parse_args([*base_args, "--combine-cis-snvs"])
    assert args_on.combine_cis_snvs is True


def test_maybe_combine_cis_snvs_opt_in_and_sample_resolution(tmp_path):
    vcf_path = tmp_path / "braf_phased.vcf"
    _write_phased_mnv_vcf(vcf_path, multi_sample=False)
    vc = load_vcf(str(vcf_path), genome="GRCh37:75")
    assert len(vc) == 2

    # Off by default: leaves SNVs unchanged.
    unchanged = maybe_combine_cis_snvs(vc, Namespace(combine_cis_snvs=False, tumor_sample_name=None))
    assert unchanged is vc
    assert len(unchanged) == 2

    # On with single-sample VCF: auto-picks sole sample and merges cis in-codon SNVs into MNV.
    combined = maybe_combine_cis_snvs(vc, Namespace(combine_cis_snvs=True, tumor_sample_name=None))
    assert len(combined) == 1
    mnv = combined[0]
    assert (mnv.contig, mnv.start, mnv.ref, mnv.alt) == ("7", 140453135, "CA", "TT")
    assert len(mnv.source_variants) == 2

    # Multi-sample VCF requires explicit --tumor-sample-name.
    multi_path = tmp_path / "braf_multi.vcf"
    _write_phased_mnv_vcf(multi_path, multi_sample=True)
    multi_vc = load_vcf(str(multi_path), genome="GRCh37:75")
    with pytest.raises(ValueError, match="--tumor-sample-name"):
        maybe_combine_cis_snvs(multi_vc, Namespace(combine_cis_snvs=True, tumor_sample_name=None))

    combined_multi = maybe_combine_cis_snvs(
        multi_vc, Namespace(combine_cis_snvs=True, tumor_sample_name="TUMOR")
    )
    assert len(combined_multi) == 1
    assert (combined_multi[0].ref, combined_multi[0].alt) == ("CA", "TT")
