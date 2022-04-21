# GROUP: global
# CONDA: cliquesnv = 2.0.3
# CONDA: samtools = 1.15.1

import subprocess
from pathlib import Path


def main(fname_bam, fname_reference, fname_marker, dname_work, seq_type, thread_num):
    dname_work.mkdir(parents=True, exist_ok=True)

    # prepare environment
    subprocess.run(
        ["samtools", "view", "-h", "-o", dname_work / "reads.sam", fname_bam]
    )

    cliquesnv_mode = None
    if seq_type == "illumina":
        cliquesnv_mode = "snv-illumina"
    elif seq_type == "pacbio":
        cliquesnv_mode = "snv-pacbio"
    else:
        raise RuntimeError(f"Invalid sequence technology: {seq_type}")

    # execute tool
    subprocess.run(
        [
            "cliquesnv",
            "-m",
            cliquesnv_mode,
            "-in",
            dname_work / "reads.sam",
            "-outDir",
            dname_work / "output",
            "-threads",
            str(thread_num),
        ]
    )

    (dname_work / "output" / "reads.fasta").rename(
        fname_marker.parent / "haplotypes.fasta"
    )


if __name__ == "__main__":
    main(
        Path(snakemake.input.fname_bam),
        Path(snakemake.input.fname_reference),
        Path(snakemake.output.fname_marker),
        Path(snakemake.output.dname_work),
        snakemake.wildcards.seq_technology,
        snakemake.threads,
    )
