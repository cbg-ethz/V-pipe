import subprocess
from pathlib import Path


def main(fname_bam, outdir):
    outdir.mkdir(parents=True, exist_ok=True)

    subprocess.run(
        [
            "tinycov",
            "covplot",
            "--skip",
            "1",
            "--ploidy",
            "0",
            "--res",
            "1",
            "--out",
            outdir / "coverage.pdf",
            fname_bam,
        ],
        check=True,
    )


if __name__ == "__main__":
    main(Path(snakemake.input.fname_bam), Path(snakemake.output.outdir))
