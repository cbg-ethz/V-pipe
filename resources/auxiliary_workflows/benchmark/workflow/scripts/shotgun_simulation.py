import re
import shutil
import tempfile
import fileinput
import subprocess
from pathlib import Path
import math

import numpy as np


RNG_SEED = 42 + int(snakemake.wildcards.replicate)
np.random.seed(RNG_SEED)


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
            "--rndSeed",
            str(RNG_SEED),
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
        assert math.isclose(
            sum(freq_list), 1
        ), f"Invalid haplotype pattern: {haplotype_pattern}, sum is {sum(freq_list)}"

    filelist_sam = []
    filelist_fastq = []

    for i, freq in enumerate(freq_list):
        haplotype_name = f"haplotype{i:04}"
        fname_haplotype = dname_work / f"{haplotype_name}.fasta"
        coverage_haplotype = int(params["coverage"] * freq)

        # simulate haplotype reads
        art_prefix = dname_work / f"art_output/haplo_{haplotype_name}_"
        art_prefix.parent.mkdir(parents=True, exist_ok=True)

        if params["seq_tech"] == "illumina":
            simulate_illumina(
                fname_haplotype, coverage_haplotype, params["read_length"], art_prefix
            )
        elif params["seq_tech"] == "pacbio":
            simulate_pacbio(
                fname_haplotype, coverage_haplotype, params["read_length"], art_prefix
            )
        elif params["seq_tech"] == "nanopore":
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
                if "@SQ" in line:
                    # replace reference name in header
                    print(
                        re.sub(
                            rf"SN:{haplotype_name}.*?\t", f"SN:{master_name}\t", line
                        ),
                        end="",
                    )
                else:
                    # only replace reference but not query name
                    print(
                        re.sub(rf"\t{haplotype_name}.*?\t", f"\t{master_name}\t", line),
                        end="",
                    )

        # gather files
        filelist_sam.append(fname_sam)
        if params["seq_tech"] == "illumina":
            reads_f = dname_work / f"art_output/haplo_{haplotype_name}_1.fq"
            reads_r = dname_work / f"art_output/haplo_{haplotype_name}_2.fq"
            filelist_fastq.append(reads_f)  # paired end reads
            filelist_fastq.append(reads_r)  # paired end reads
        elif params["seq_tech"] in ["pacbio", "nanopore"]:
            reads_f = dname_work / f"art_output/haplo_{haplotype_name}__0001.fastq"
            filelist_fastq.append(reads_f)

    # merge sam files
    fname_merged_sam = dname_work / "haplotype_reads.sam"
    subprocess.run(["samtools", "merge", "-f", "-o", fname_merged_sam, *filelist_sam])

    # save result
    with open(fname_fastq, "wb") as fd_write:
        for fname in filelist_fastq:
            with open(fname, "rb") as fd_read:
                shutil.copyfileobj(fd_read, fd_write)

    with tempfile.NamedTemporaryFile() as fd_tmp:
        subprocess.run(
            ["samtools", "view", "-b", "-h", "-o", fd_tmp.name, fname_merged_sam]
        )
        subprocess.run(["samtools", "sort", "-o", fname_bam, fd_tmp.name])


if __name__ == "__main__":
    main(
        Path(snakemake.output.fname_fastq),
        Path(snakemake.output.fname_bam),
        Path(snakemake.input.dname_work),
        snakemake.params.haplotype_generation,
        snakemake.params.params,
    )
