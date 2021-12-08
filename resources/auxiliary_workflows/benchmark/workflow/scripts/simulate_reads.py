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
    # initial setup
    np.random.seed(42)
    dname_work.mkdir(parents=True, exist_ok=True)

    # generate random master sequence
    seq_master = "".join(np.random.choice(BASE_LIST, size=1000))
    fname_reference.write_text(f">MasterSequence\n{seq_master}\n")

    # infer haplotype sequences
    haplotype_count = 10
    mutation_rate = 0.2
    insertion_rate = 0.1
    deletion_rate = 0.05

    haplotype_dict = {}
    for i in range(haplotype_count):
        seq_haplotype = generate_haplotype(
            seq_master, mutation_rate, insertion_rate, deletion_rate
        )
        haplotype_dict[f"Haplotype_{i:04}"] = seq_haplotype

    fname_haplotypes = dname_work / "haplotypes.fasta"
    fname_haplotypes.write_text(
        "\n".join(f">{k}\n{v}" for k, v in haplotype_dict.items())
    )

    # simulate reads for each haplotype
    art_prefix = dname_work / "art_output"
    subprocess.run(
        [
            "art_illumina",
            "-sam",
            "-i",
            fname_haplotypes,
            "-c",
            "100",
            "-l",
            "250",
            "-o",
            art_prefix,
        ]
    )

    # save result
    art_prefix.with_suffix(".fq").rename(fname_fastq)
    subprocess.run(
        ["samtools", "view", "-b", "-o", fname_bam, art_prefix.with_suffix(".sam")]
    )


if __name__ == "__main__":
    main(
        Path(snakemake.output.fname_fastq),
        Path(snakemake.output.fname_bam),
        Path(snakemake.output.fname_reference),
        Path(snakemake.output.dname_work),
    )
