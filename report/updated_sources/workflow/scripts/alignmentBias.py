#!/usr/bin/env python3

import argparse
import pysam
from Bio import SeqIO, pairwise2

__author__ = "Susana Posada-Cespedes"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"


def parse_args():
    """Set up the parsing of command-line arguments"""

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    requiredNamed = parser.add_argument_group("required named arguments")
    requiredNamed.add_argument(
        "-r",
        required=True,
        default=None,
        metavar="FASTA",
        dest="reference_file",
        type=str,
        help="Referene sequence",
    )
    requiredNamed.add_argument(
        "-b",
        required=True,
        default=None,
        metavar="BAM",
        dest="bamfile",
        help="BAM file",
    )
    requiredNamed.add_argument(
        "-f",
        required=True,
        default=None,
        metavar="FASTQ",
        dest="fastq_file",
        help="FASTQ file containing reads after QC",
    )
    requiredNamed.add_argument(
        "--hap",
        required=True,
        default=None,
        metavar="FASTA",
        dest="haplotype_seqs",
        help="FASTA file containing haplotype sequences",
    )
    requiredNamed.add_argument(
        "-N",
        required=False,
        default="sample",
        metavar="STR",
        dest="sampleID",
        help="Patient/sample identifiers",
    )
    parser.add_argument(
        "-p",
        required=False,
        default=False,
        action="store_true",
        dest="paired",
        help="Indicate whether to simulate paired-end reads",
    )
    parser.add_argument(
        "-o",
        required=False,
        default="alignment_bias.tsv",
        metavar="FILENAME",
        dest="output_file",
        type=str,
        help="Output file",
    )

    return parser.parse_args()


def hamming_dist(s1, s2):
    """Return the Hamming distance between sequences"""
    assert len(s1) == len(s2), (
        "Hamming distance is undefined for sequences " "differing in their lengths"
    )
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))


def main():
    args = parse_args()
    reference = SeqIO.read(args.reference_file, "fasta")
    reference_len = len(reference)
    start = None
    end = None

    hap_cnt = {}
    hap_coverage = {}
    with pysam.AlignmentFile(args.bamfile, "rb") as alnfile:
        for read in alnfile.fetch(reference=reference.id, start=start, end=end):
            query_name = read.query_name
            hapID = query_name.split("-")[0]
            if hapID in hap_cnt:
                hap_cnt[hapID] += 1
            else:
                hap_cnt[hapID] = 1

            # Get the percent of aligned bases per read
            cigar_arr = read.get_cigar_stats()[0]
            aux = cigar_arr[0] + cigar_arr[1]
            assert aux == read.query_alignment_length
            percent_aligned = read.query_alignment_length / read.infer_read_length()
            assert percent_aligned <= 1, f"{cigar_arr}"
            if hapID in hap_coverage:
                hap_coverage[hapID] += percent_aligned
            else:
                hap_coverage[hapID] = percent_aligned

    for k, v in hap_cnt.items():
        aux = hap_coverage[k] / v
        assert aux <= 1, f"{aux}"
        hap_coverage[k] = aux

    # Parse FASTQ file to compute the expected number of reads per haplotype
    hap_expected = {}
    total = 0
    for read in SeqIO.parse(args.fastq_file, "fastq"):
        total += 1
        hapID = read.id.split("-")[0]
        if hapID in hap_expected:
            hap_expected[hapID] += 1
        else:
            hap_expected[hapID] = 1

    percent_aligned = {}
    for k, v in hap_expected.items():
        if k in hap_cnt:
            if args.paired:
                percent_aligned[k] = hap_cnt[k] / (v * 2)
            else:
                percent_aligned[k] = hap_cnt[k] / v
        else:
            percent_aligned[k] = 0
            hap_coverage[k] = 0

    # Compute divergence from the reference sequence
    divergence = {}
    affine_pen = pairwise2.affine_penalty(-1, -0.1, True)
    with open(args.haplotype_seqs, "r") as infile:
        for record in SeqIO.parse(infile, "fasta"):
            alignment = pairwise2.align.globalmc(
                reference.seq,
                record.seq,
                1,
                0,
                affine_pen,
                affine_pen,
                one_alignment_only=True,
            )
            divergence[record.id] = (
                hamming_dist(alignment[0][0], alignment[0][1]) / alignment[0][4]
            )

    with open(args.output_file, "w") as outfile:
        outfile.write(
            "SampleID\tHaplotypeID\tDivergence\tPercent-aligned\t"
            "Percent-bases-aligned\n"
        )
        for k, v in percent_aligned.items():
            outfile.write(
                f"{args.sampleID}\t{k}\t{divergence[k]}\t{v}\t" f"{hap_coverage[k]}\n"
            )


if __name__ == "__main__":
    main()
