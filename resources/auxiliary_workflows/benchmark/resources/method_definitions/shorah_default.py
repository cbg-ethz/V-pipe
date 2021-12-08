import subprocess
from pathlib import Path


def main(fname_bam, fname_reference, fname_results, dname_work):
    # TODO: implement shorah
    # subprocess.run(
    #     ["shorah", "snv", "-b", fname_bam, "-f", fname_reference, "--out_format=vcf"]
    # )

    dname_work.mkdir(parents=True, exist_ok=True)
    fname_results.touch()


if __name__ == "__main__":
    main(
        Path(snakemake.input.fname_bam),
        Path(snakemake.input.fname_reference),
        Path(snakemake.output.fname_results),
        Path(snakemake.output.dname_work),
    )
