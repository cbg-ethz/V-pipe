# GROUP: global
# CONDA: quasirecomb = 1.2

import subprocess
from pathlib import Path


def main(fname_bam, fname_reference, fname_marker, dname_work):
    dname_work.mkdir(parents=True, exist_ok=True)

    subprocess.run(["samtools", "index", fname_bam])
    subprocess.run(["quasirecomb", "-i", fname_bam, "-o", dname_work / "output"])

    (dname_work / "output" / "quasispecies.fasta").rename(
        fname_marker.parent / "haplotypes.fasta"
    )


if __name__ == "__main__":
    main(
        Path(snakemake.input.fname_bam),
        Path(snakemake.input.fname_reference),
        Path(snakemake.output.fname_marker),
        Path(snakemake.output.dname_work),
    )
