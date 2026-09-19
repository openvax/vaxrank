"""Account for every PacBio groupdedup molecule before/after published mapping."""

import argparse
from collections import Counter
from pathlib import Path
import json
import subprocess

import pysam

from build import digest, now, write_json


def molecule_number(read):
    prefix, number = read.query_name.split("/")
    if prefix != "molecule" or not number.isdigit():
        raise ValueError("Unexpected PacBio query-name scheme: " + read.query_name)
    return int(number)


def main(root):
    directory = root / "pacbio_full"
    mapped = directory / "IPISRC044_T1_sclrs_live_pbmm2_mapped.bam"
    original = directory / "IPISRC044_T1_sclrs_live.seg.lima.flt.ref.corr.sort.dedup.bam"
    mapped_names = bytearray()
    counts = Counter()
    length_histograms = {k: Counter() for k in ("input", "mapped_input", "absent_from_mapped")}
    with pysam.AlignmentFile(mapped) as bam:
        programs = bam.header.to_dict().get("PG", [])
        for read in bam:
            number = molecule_number(read)
            if number >= len(mapped_names):
                mapped_names.extend(bytes(number + 1 - len(mapped_names)))
            counts["mapped_bam_records"] += 1
            counts["unique_mapped_names"] += not mapped_names[number]
            mapped_names[number] = 1
            counts["secondary_records"] += read.is_secondary
            counts["supplementary_records"] += read.is_supplementary
            counts["low_mapq_records"] += read.mapping_quality < 20
    print("Mapped inventory:", dict(counts), flush=True)
    seen_input = bytearray(len(mapped_names))
    with pysam.AlignmentFile(original, check_sq=False) as bam, pysam.AlignmentFile(
            directory / "absent-from-published-mapping.bam", "wb", header=bam.header) as missing, (
            directory / "absent-from-published-mapping.fa").open("w") as fasta:
        for read in bam:
            number = molecule_number(read)
            if number >= len(seen_input):
                seen_input.extend(bytes(number + 1 - len(seen_input)))
            if seen_input[number]:
                raise ValueError("Nonunique deduplicated molecule name")
            seen_input[number] = 1
            counts["input_molecules"] += 1
            length_histograms["input"][read.query_length] += 1
            present = number < len(mapped_names) and mapped_names[number]
            length_histograms["mapped_input" if present else "absent_from_mapped"][read.query_length] += 1
            if present:
                counts["input_molecules_represented_in_mapping"] += 1
            else:
                counts["input_molecules_absent_from_mapping"] += 1
                missing.write(read)
                if read.query_sequence:
                    fasta.write(">%s\n%s\n" % (read.query_name, read.get_forward_sequence()))
        if any(flag and not seen_input[i] for i, flag in enumerate(mapped_names)):
            raise ValueError("Mapped BAM contains molecule identities absent from input")
    output = dict(created_utc=now(), counts=dict(counts), programs=programs,
                  original_sha256=digest(original), mapped_sha256=digest(mapped),
                  absent_bam_sha256=digest(directory / "absent-from-published-mapping.bam"),
                  absent_fasta_sha256=digest(directory / "absent-from-published-mapping.fa"),
                  length_histograms=length_histograms,
                  interpretation="Identity accounting only: missing molecules may fail alignment/quality criteria; any remapping requires independent evidence gates")
    write_json(directory / "molecule-accounting.json", output)
    print(dict(counts), flush=True)


def realign(root, reference):
    directory = root / "pacbio_full"
    reads = directory / "absent-from-published-mapping.fa"
    output = directory / "rescued-genome.bam"
    command = ["minimap2", "-ax", "splice:hq", "-Y", "-t", "4", str(reference), str(reads)]
    sort_command = ["samtools", "sort", "-@", "2", "-o", str(directory / "rescued-before-tags.bam"), "-"]
    with (directory / "realignment.stderr.log").open("w") as log:
        alignment = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=log)
        sorting = subprocess.Popen(sort_command, stdin=alignment.stdout, stderr=log)
        alignment.stdout.close()
        sort_status = sorting.wait()
        align_status = alignment.wait()
        if align_status or sort_status:
            raise RuntimeError("Alignment/sort failed: %s/%s" % (align_status, sort_status))
    with pysam.AlignmentFile(directory / "absent-from-published-mapping.bam", check_sq=False) as bam:
        originals = {r.query_name: r for r in bam}
        groups = bam.header.to_dict().get("RG", [])
    counts = Counter()
    with pysam.AlignmentFile(directory / "rescued-before-tags.bam") as bam:
        header = bam.header.to_dict()
        header["RG"] = groups
        with pysam.AlignmentFile(output, "wb", header=header) as target:
            for read in bam:
                original = originals[read.query_name]
                if read.get_forward_sequence() != original.get_forward_sequence():
                    raise ValueError("Remapping changed a stored molecule sequence")
                for tag, value, value_type in original.get_tags(with_value_type=True):
                    if value_type == "B":
                        read.set_tag(tag, value)
                    else:
                        read.set_tag(tag, value, value_type=value_type)
                qualities = original.get_forward_qualities()
                read.query_qualities = qualities[::-1] if qualities is not None and read.is_reverse else qualities
                target.write(read)
                counts["records"] += 1
                counts["mapped_primary"] += not (read.is_unmapped or read.is_secondary or read.is_supplementary)
                counts["unmapped"] += read.is_unmapped
    pysam.index(str(output))
    write_json(directory / "realignment.json", dict(command=command, sort_command=sort_command,
        reference=str(reference), reference_sha256=digest(reference), input_sha256=digest(reads),
        output_sha256=digest(output), counts=dict(counts), minimap2=subprocess.check_output(["minimap2", "--version"], text=True).strip(),
        source_molecule_provenance="Original CB/XM/im/is/ic and other input tags restored by exact unique molecule identity; no independent library added"))
    print("Remapped omitted molecules:", json.dumps(counts), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--realign", action="store_true")
    parser.add_argument("--reference", type=Path)
    args = parser.parse_args()
    if args.realign:
        if args.reference is None:
            parser.error("--realign requires --reference")
        realign(args.run, args.reference)
    else:
        main(args.run)
