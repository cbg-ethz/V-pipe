# GROUP: local
# CONDA: boost = 1.77.0
# CONDA: htslib = 1.14
# PIP: git+https://github.com/LaraFuhrmann/shorah@master

import subprocess
from pathlib import Path


def main(fname_bam, fname_reference, fname_insert_bed, fname_marker, dname_work):
    dname_work.mkdir(parents=True, exist_ok=True)

    if fname_insert_bed == "":
        # no insert file --> shotgun mode
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
    else:
        # insert file --> amplicon mode
        subprocess.run(
            [
                "shorah",
                "shotgun",
                "-b",
                fname_bam.resolve(),
                "-f",
                fname_reference.resolve(),
                "--insert-file",
                fname_insert_bed.resolve(),
            ],
            cwd=dname_work,
            check=True,
        )

    (dname_work / "snv" / "SNVs_0.010000_final.vcf").rename(
        fname_marker.parent / "snvs.vcf"
    )


if __name__ == "__main__":
    main(
        Path(snakemake.input.fname_bam),
        Path(snakemake.input.fname_reference),
        Path(snakemake.input.fname_insert_bed),
        Path(snakemake.output.fname_marker),
        Path(snakemake.output.dname_work),
    )
