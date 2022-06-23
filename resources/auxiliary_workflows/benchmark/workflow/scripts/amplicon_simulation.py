import re
import shutil
import tempfile
import fileinput
import subprocess
from pathlib import Path

import numpy as np

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import math


def cut_amplicon_regions(fname_reference, fname_insert_bed, fname_output):
    amplicons = []
    # open insert.bed file
    with open(fname_insert_bed) as f:
        for line in f:
            L = line.strip().split()
            # 1-based, exclusive
            amplicons.append((int(L[1]), int(L[2])))

    fasta_sequences = SeqIO.parse(open(fname_reference), "fasta")
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)

    records = []
    for amplicon in amplicons:
        records.append(
            SeqRecord(
                Seq(sequence[amplicon[0] - 1 : amplicon[1] - 1]),
                id=f"{name}_position{amplicon[0]}..{amplicon[1]}",
                description="",
            )
        )
    SeqIO.write(records, fname_output, "fasta")


def simulate_illumina(fname_haplotype, coverage_haplotype, read_length, art_prefix):
    print(fname_haplotype)
    subprocess.run(
        [
            "art_illumina",
            "-sam",
            "--paired",
            "-m",
            str(
                read_length * 2 * 0.8
            ),  # mean and standard deviation of DNA/RNA fragment lengths
            "-s",
            str(10),
            "-i",
            fname_haplotype,
            "-f",
            str(coverage_haplotype),
            "-l",
            str(read_length),
            "-o",
            art_prefix,
        ]
    )


def main(
    fname_insert_bed,
    dname_work,
    fname_fastq_R1,
    fname_fastq_R2,
    fname_bam,
    haplotype_generation,
    params,
):
    master_name = "MasterSequence"

    if haplotype_generation == "distance":
        (
            n_group1,
            n_group2,
            d_group12,
            d_group1,
            d_group2,
            freq_dist,
            freq_param,
        ) = params["haplos"].split("@")
        n_group1 = int(n_group1)
        n_group2 = int(n_group2)
        n_haplo = n_group1 + n_group2
        # obtain haplotype frequencies
        if freq_dist == "dirichlet":
            if type(freq_param) == str:
                alpha = [float(freq) for freq in freq_param.split(":")]
            if type(freq_param) == np.float64:
                if np.isnan(freq_param):
                    alpha = np.ones(n_haplo)

            freq_list = np.random.dirichlet(alpha)

        elif freq_dist == "geom":
            freq_param = float(freq_param)
            freq_list = np.asarray([freq_param ** (i + 1) for i in range(n_haplo)])
            freq_list = freq_list / np.sum(freq_list)

    # infer haplotype sequences
    elif haplotype_generation == "mutation_rate":
        haplotype_pattern = params["haplos"].split("@")[-1]
        # infer haplotype sequences
        freq_list = [float(freq) for freq in haplotype_pattern.split(":")]
        assert (
            math.isclose(sum(freq_list),1)
        ), f"Invalid haplotype pattern: {haplotype_pattern}, sum is {sum(freq_list)}"

    filelist_sam = []
    filelist_fastq_f = []
    filelist_fastq_r = []

    for i, freq in enumerate(freq_list):
        haplotype_name = f"haplotype{i:04}"
        fname_haplotype = dname_work / f"{haplotype_name}.fasta"
        fname_insert_haplotype = dname_work / f"{haplotype_name}.insert.fasta"
        cut_amplicon_regions(fname_haplotype, fname_insert_bed, fname_insert_haplotype)

        coverage_haplotype = int(params["coverage"] * freq)

        # simulate haplotype reads
        art_prefix = dname_work / f"art_output/haplo_{haplotype_name}_"
        art_prefix.parent.mkdir(parents=True, exist_ok=True)

        if params["seq_tech"] == "illumina":
            simulate_illumina(
                fname_insert_haplotype,
                coverage_haplotype,
                params["read_length"],
                art_prefix,
            )

        fname_sam = art_prefix.with_suffix(".sam")
        with fileinput.FileInput(fname_sam, inplace=True, backup=".bak") as fd:
            need_header = True
            for line in fd:
                if line.startswith("@SQ"):
                    # only keep a single @SQ header line
                    if need_header:
                        print(
                            f"@SQ\tSN:{master_name}\tLN:{snakemake.wildcards.genome_size}"
                        )
                        need_header = False
                else:
                    print(
                        re.sub(rf"{haplotype_name}_.*?\t", f"{master_name}\t", line),
                        end="",
                    )

        # gather files
        filelist_sam.append(fname_sam)
        if params["seq_tech"] == "illumina":
            reads_f = dname_work / f"art_output/haplo_{haplotype_name}_1.fq"
            reads_r = dname_work / f"art_output/haplo_{haplotype_name}_2.fq"
            filelist_fastq_f.append(reads_f)  # paired end reads
            filelist_fastq_r.append(reads_r)  # paired end reads

    # save result
    with open(fname_fastq_R1, "wb") as fd_write:
        for fname in filelist_fastq_f:
            with open(fname, "rb") as fd_read:
                shutil.copyfileobj(fd_read, fd_write)

    # save result
    with open(fname_fastq_R2, "wb") as fd_write:
        for fname in filelist_fastq_r:
            with open(fname, "rb") as fd_read:
                shutil.copyfileobj(fd_read, fd_write)

    # convert merged sam files to bam
    fname_merged_sam = dname_work / "haplotype_reads.sam"
    subprocess.run(["samtools", "merge", "-f", "-o", fname_merged_sam, *filelist_sam])

    with tempfile.NamedTemporaryFile() as fd_tmp:
        subprocess.run(
            ["samtools", "view", "-b", "-h", "-o", fd_tmp.name, fname_merged_sam]
        )
        subprocess.run(["samtools", "sort", "-o", fname_bam, fd_tmp.name])


if __name__ == "__main__":
    main(
        Path(snakemake.input.fname_insert_bed),
        Path(snakemake.input.dname_work),
        Path(snakemake.output.fname_fastq_R1),
        Path(snakemake.output.fname_fastq_R2),
        Path(snakemake.output.fname_bam),
        snakemake.params.haplotype_generation,
        snakemake.params.params,
    )
