# GROUP: local
# CONDA: shorah = 1.99.2

import subprocess
from pathlib import Path


def main(fname_bam, fname_reference, fname_result, fname_result_haplos, dname_work, seq_tech):
    dname_work.mkdir(parents=True, exist_ok=True)

    if seq_tech == 'illumina':
        window_length = 162

        subprocess.run(
            [
                "shorah",
                "shotgun",
                "-b",
                fname_bam.resolve(),
                "-f",
                fname_reference.resolve(),
                "--windowsize",
                str(window_length),
            ],
            cwd=dname_work,
            check=True,
        )
    else:
        subprocess.run(
            [
                "shorah",
                "shotgun",
                "-b",
                fname_bam.resolve(),
                "-f",
                fname_reference.resolve(),
            ],
            cwd=dname_work,
            check=True,
        )
    (dname_work / "snv" / "SNVs_0.010000_final.vcf").rename(fname_result)

    open(fname_result_haplos, "a").close()


if __name__ == "__main__":
    main(
        Path(snakemake.input.fname_bam),
        Path(snakemake.input.fname_reference),
        Path(snakemake.output.fname_result),
        Path(snakemake.output.fname_result_haplos),
        Path(snakemake.output.dname_work),
        snakemake.wildcards.seq_tech,
    )
