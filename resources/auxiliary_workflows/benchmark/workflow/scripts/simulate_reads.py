import shutil
import subprocess
from pathlib import Path

import numpy as np


BASE_LIST = list("TCGA")


def generate_haplotype(seq_master, mutation_rate=0, insertion_rate=0, deletion_rate=0):
    """Generate haplotype from master sequence."""
    seq_haplotype = np.asarray(list(seq_master))

    # deletions
    deletion_count = int(len(seq_haplotype) * deletion_rate)
    position_list = np.random.choice(
        np.arange(len(seq_haplotype)), size=deletion_count, replace=False
    )
    seq_haplotype = np.delete(seq_haplotype, position_list)

    # mutations
    mutation_count = int(len(seq_haplotype) * mutation_rate)
    position_list = np.random.choice(
        np.arange(len(seq_haplotype)), size=mutation_count, replace=False
    )
    seq_haplotype[position_list] = np.random.choice(
        BASE_LIST, size=len(position_list)
    )  # TODO: ensure that mutated base is not the same as original base

    # insertions
    insertion_count = int(len(seq_haplotype) * insertion_rate)
    position_list = np.random.choice(
        np.arange(len(seq_haplotype)), size=insertion_count, replace=False
    )
    seq_haplotype = np.insert(
        seq_haplotype,
        position_list,
        np.random.choice(BASE_LIST, size=insertion_count),
    )

    return "".join(seq_haplotype)


def main(fname_fastq, fname_bam, fname_reference, dname_work):
    """Create master sequence, infer haplotypes and simulate reads."""
    genome_size = 1000
    coverage = 100
    read_length = 250
    haplotype_pattern = "0.5;0.2;0.3"
    mutation_rate = 0.2
    insertion_rate = 0.1
    deletion_rate = 0.05

    # initial setup
    np.random.seed(42)
    dname_work.mkdir(parents=True, exist_ok=True)

    # generate random master sequence
    seq_master = "".join(np.random.choice(BASE_LIST, size=genome_size))
    fname_reference.write_text(f">MasterSequence\n{seq_master}\n")

    # infer haplotype sequences
    freq_list = [float(freq) for freq in haplotype_pattern.split(";")]
    assert sum(freq_list) == 1, f"Invalid haplotype pattern: {haplotype_pattern}"

    filelist_sam = []
    filelist_fastq = []
    for i, freq in enumerate(freq_list):
        # generate haplotype
        haplotype_name = f"haplotype_{i:04}"
        seq_haplotype = generate_haplotype(
            seq_master, mutation_rate, insertion_rate, deletion_rate
        )

        coverage_haplotype = int(coverage * freq)

        # save haplotype in FASTA
        fname_haplotype = dname_work / f"{haplotype_name}.fasta"
        fname_haplotype.write_text(f">{haplotype_name}\n{seq_haplotype}\n")

        # simulate haplotype reads
        art_prefix = dname_work / f"art_output/haplo_{haplotype_name}"
        art_prefix.parent.mkdir(parents=True, exist_ok=True)

        subprocess.run(
            [
                "art_illumina",
                "-sam",
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

        # gather files
        filelist_sam.append(art_prefix.with_suffix(".sam"))
        filelist_fastq.append(art_prefix.with_suffix(".fq"))

    # merge sam files
    fname_merged_sam = dname_work / "haplotype_reads.sam"
    subprocess.run(["samtools", "merge", "-o", fname_merged_sam, *filelist_sam])

    # save result
    with open(fname_fastq, "wb") as fd_write:
        for fname in filelist_fastq:
            with open(fname, "rb") as fd_read:
                shutil.copyfileobj(fd_read, fd_write)

    subprocess.run(["samtools", "view", "-b", "-o", fname_bam, fname_merged_sam])


if __name__ == "__main__":
    main(
        Path(snakemake.output.fname_fastq),
        Path(snakemake.output.fname_bam),
        Path(snakemake.output.fname_reference),
        Path(snakemake.output.dname_work),
    )
