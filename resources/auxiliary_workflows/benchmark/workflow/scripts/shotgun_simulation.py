import shutil
import tempfile
import fileinput
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd


def simulate_illumina(fname_haplotype, coverage_haplotype, read_length, art_prefix):
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
            "-c",
            str(coverage_haplotype),
            "-l",
            str(read_length),
            "-o",
            art_prefix,
        ]
    )


def simulate_pacbio(fname_haplotype, coverage_haplotype, read_length, art_prefix):
    """Simulate long reads with PacBio error profile using PBSIM2."""
    diff_ratio = "6:50:54"  # --difference-ratio
    hmm_model = "resources/pbsim2/data/P6C4.model"

    # run PBSIM2
    subprocess.run(
        [
            "pbsim",
            fname_haplotype,
            "--length-min",
            str(read_length - 5),
            "--length-max",
            str(read_length),
            "--length-mean",
            str(read_length - 2),
            "--depth",
            str(coverage_haplotype),
            "--hmm_model",
            hmm_model,
            "--difference-ratio",
            diff_ratio,
            "--prefix",
            art_prefix,
        ]
    )
    # align reads
    subprocess.run(["bwa", "index", fname_haplotype])

    subprocess.run(
        [
            "bwa",
            "mem",
            "-o",
            str(art_prefix) + ".sam",
            fname_haplotype,
            str(art_prefix) + "_0001.fastq",
        ]
    )

    subprocess.run(
        ["samtools", "sort", str(art_prefix) + ".sam", "-o", str(art_prefix) + ".sam"]
    )


def simulate_nanopore(fname_haplotype, coverage_haplotype, read_length, art_prefix):
    # --difference-ratio 23:31:46
    diff_ratio = "23:31:46"  # --difference-ratio
    hmm_model = "resources/pbsim2/data/P6C4.model"

    # run PBSIM2
    subprocess.run(
        [
            "pbsim",
            fname_haplotype,
            "--length-min",
            str(read_length - 5),
            "--length-max",
            str(read_length),
            "--length-mean",
            str(read_length - 2),
            "--depth",
            str(coverage_haplotype),
            "--hmm_model",
            hmm_model,
            "--difference-ratio",
            diff_ratio,
            "--prefix",
            art_prefix,
        ]
    )
    # align reads
    subprocess.run(["bwa", "index", fname_haplotype])

    subprocess.run(
        [
            "bwa",
            "mem",
            "-o",
            str(art_prefix) + ".sam",
            fname_haplotype,
            str(art_prefix) + "_0001.fastq",
        ]
    )

    subprocess.run(
        ["samtools", "sort", str(art_prefix) + ".sam", "-o", str(art_prefix) + ".sam"]
    )


def main(fname_fastq, fname_bam, dname_work, haplotype_generation, params):

    master_name = "MasterSequence"

    if haplotype_generation == "distance":
        n_haplo = params["n_group1"] + params["n_group2"]
        # obtain haplotype frequencies
        if params["freq_distribution"] == "dirichlet":
            if type(params["freq_param"]) == str:
                alpha = [float(freq) for freq in params["freq_param"].split(":")]
            if type(params["freq_param"]) == np.float64:
                if np.isnan(params["freq_param"]):
                    alpha = np.ones(n_haplo)

            freq_list = np.random.dirichlet(alpha)

        elif params["freq_distribution"] == "geom":
            freq_param = float(params["freq_param"])
            freq_list = np.asarray([freq_param ** (i + 1) for i in range(n_haplo)])
            freq_list = freq_list / np.sum(freq_list)

    # infer haplotype sequences
    elif haplotype_generation == "mutation_rate":
        # infer haplotype sequences
        freq_list = [float(freq) for freq in params["haplotype_pattern"].split(":")]
        assert (
            sum(freq_list) == 1
        ), f"Invalid haplotype pattern: {params['haplotype_pattern']}, sum is {sum(freq_list)}"

    filelist_sam = []
    filelist_fastq = []

    for i, freq in enumerate(freq_list):
        haplotype_name = f"haplotype{i:04}"
        fname_haplotype = dname_work / f"{haplotype_name}.fasta"
        coverage_haplotype = int(params["coverage"] * freq)

        # simulate haplotype reads
        art_prefix = dname_work / f"art_output/haplo_{haplotype_name}_"
        art_prefix.parent.mkdir(parents=True, exist_ok=True)

        if params["seq_technology"] == "illumina":
            simulate_illumina(
                fname_haplotype, coverage_haplotype, params["read_length"], art_prefix
            )
        elif params["seq_technology"] == "pacbio":
            simulate_pacbio(
                fname_haplotype, coverage_haplotype, params["read_length"], art_prefix
            )
        elif params["seq_technology"] == "nanopore":
            simulate_nanopore(
                fname_haplotype, coverage_haplotype, params["read_length"], art_prefix
            )

        # we assume that generated SAM file represent reads being mapped
        # to a reference (or de-novo consensus sequence).
        # We must thus change the reference name accordingly
        # TODO: there is probably/hopefully a better way of doing this
        fname_sam = art_prefix.with_suffix(".sam")
        with fileinput.FileInput(fname_sam, inplace=True, backup=".bak") as fd:
            for line in fd:
                print(line.replace(haplotype_name, master_name), end="")

        # gather files
        filelist_sam.append(fname_sam)
        if params["seq_technology"] == "illumina":
            reads_f = dname_work / f"art_output/haplo_{haplotype_name}_1.fq"
            reads_r = dname_work / f"art_output/haplo_{haplotype_name}_2.fq"
            filelist_fastq.append(reads_f)  # paired end reads
            filelist_fastq.append(reads_r)  # paired end reads
        elif params["seq_technology"] in ["pacbio", "nanopore"]:
            reads_f = dname_work / f"art_output/haplo_{haplotype_name}__0001.fastq"
            filelist_fastq.append(reads_f)

    # merge sam files
    fname_merged_sam = dname_work / "haplotype_reads.sam"
    subprocess.run(["samtools", "merge", "-o", fname_merged_sam, *filelist_sam])

    # save result
    with open(fname_fastq, "wb") as fd_write:
        for fname in filelist_fastq:
            with open(fname, "rb") as fd_read:
                shutil.copyfileobj(fd_read, fd_write)

    with tempfile.NamedTemporaryFile() as fd_tmp:
        subprocess.run(["samtools", "view", "-b", "-o", fd_tmp.name, fname_merged_sam])
        subprocess.run(["samtools", "sort", "-o", fname_bam, fd_tmp.name])


if __name__ == "__main__":
    main(
        Path(snakemake.output.fname_fastq),
        Path(snakemake.output.fname_bam),
        Path(snakemake.input.dname_work),
        snakemake.params.haplotype_generation,
        snakemake.params.params,
    )
