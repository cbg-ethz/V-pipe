# GROUP: local
# CONDA: lofreq = 2.1.5

import subprocess
from pathlib import Path


def main(fname_bam, fname_reference, fname_marker, dname_work):
    dname_work.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "lofreq",
            "call",
            "-f",
            fname_reference.resolve(),
            "-o",
            (fname_marker.parent / "snvs.vcf").resolve(),
            fname_bam.resolve(),
        ],
        cwd=dname_work,
    )


if __name__ == "__main__":
    main(
        Path(snakemake.input.fname_bam),
        Path(snakemake.input.fname_reference),
        Path(snakemake.output.fname_marker),
        Path(snakemake.output.dname_work),
    )
