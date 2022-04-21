# GROUP: global
# CONDA: pehaplo = v0.1_dev05
# CONDA: samtools = 1.15.1

import subprocess
from pathlib import Path


def main(fname_bam, fname_reference, fname_marker, dname_work):
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
            "pehaplo.py",
            "-f1",
            (dname_work / "reads.R1.fastq").resolve(),
            "-f2",
            (dname_work / "reads.R2.fastq").resolve(),
            "--overlap_len",
            str(int(snakemake.wildcards.read_length) // 4),
            "--read_len",
            snakemake.wildcards.read_length,
        ],
        cwd=dname_work,
        check=True,
    )

    # TODO: select correct result
    print(list(dname_work.iterdir()))


if __name__ == "__main__":
    main(
        Path(snakemake.input.fname_bam),
        Path(snakemake.input.fname_reference),
        Path(snakemake.output.fname_marker),
        Path(snakemake.output.dname_work),
    )
