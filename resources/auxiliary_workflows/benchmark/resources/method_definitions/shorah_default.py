# GROUP: local
# CONDA: shorah = 1.99.2

import subprocess
from pathlib import Path


def main(fname_bam, fname_reference, fname_marker, dname_work):
    dname_work.mkdir(parents=True, exist_ok=True)
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
    )

    (dname_work / "snv" / "SNVs_0.010000_final.vcf").rename(
        fname_marker.parent / "snvs.vcf"
    )


if __name__ == "__main__":
    main(
        Path(snakemake.input.fname_bam),
        Path(snakemake.input.fname_reference),
        Path(snakemake.output.fname_marker),
        Path(snakemake.output.dname_work),
    )
