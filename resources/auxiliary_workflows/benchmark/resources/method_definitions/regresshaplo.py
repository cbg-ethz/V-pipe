# GROUP: global
# CONDA: r-regresshaplo = 1.0

# installation crashes as dependencies are not available

import subprocess
from pathlib import Path


def main(fname_bam, fname_reference, fname_result, dname_work):
    dname_work.mkdir(parents=True, exist_ok=True)

    subprocess.run(["quasirecomb", "-i", fname_bam, "-o", dname_work / "output"])

    # TODO: aggregate results
    # (dname_work / "output" / "quasispecies.fasta").rename(fname_result)


if __name__ == "__main__":
    main(
        Path(snakemake.input.fname_bam),
        Path(snakemake.input.fname_reference),
        Path(snakemake.output.fname_result),
        Path(snakemake.output.dname_work),
    )
