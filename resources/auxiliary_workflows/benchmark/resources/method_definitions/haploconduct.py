# GROUP: global
# CONDA: haploconduct = 0.2.1
# CONDA: samtools = 1.15.1

import subprocess
from pathlib import Path


def main(fname_bam, fname_reference, fname_result, dname_work):
    dname_work.mkdir(parents=True, exist_ok=True)

    subprocess.run(
        [
            "samtools",
            "fastq",
            "-1",
            dname_work / "reads.R1.fastq",
            "-2",
            dname_work / "reads.R2.fastq",
            fname_bam,
        ],
        check=True,
    )
    subprocess.run(
        [
            "haploconduct",
            "savage",
            # "--ref",
            # fname_reference.resolve(),
            "--split",
            "4",  # TODO: choose this appropriately
            "-p1",
            (dname_work / "reads.R1.fastq").resolve(),
            "-p2",
            (dname_work / "reads.R2.fastq").resolve(),
        ],
        cwd=dname_work,
    )

    # select result (later stage is better)
    result_a = dname_work / "contigs_stage_a.fasta"
    result_b = dname_work / "contigs_stage_b.fasta"
    result_c = dname_work / "contigs_stage_c.fasta"

    for result in [result_c, result_b, result_a]:
        if result.exists():
            print(f"Using stage result {result}")
            result.rename(fname_result)
            break
    else:
        raise RuntimeError("Savage crashed and no results were generated")


if __name__ == "__main__":
    main(
        Path(snakemake.input.fname_bam),
        Path(snakemake.input.fname_reference),
        Path(snakemake.output.fname_result),
        Path(snakemake.output.dname_work),
    )
