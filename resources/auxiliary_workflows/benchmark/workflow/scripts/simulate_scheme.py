import subprocess
from pathlib import Path

def main(fname_reference, fname_insert_bed, dname_out, params):

    amplicon_size = params['seq_mode'].split(":")[1]
    target_overlap = params['seq_mode'].split(":")[2]

    subprocess.run(
        [
            "primalscheme",
            "multiplex",
            "--amplicon-size",
            str(amplicon_size),
            "--name",
            "reference",
            "--target-overlap",
            str(target_overlap),
            "--force",
            "--outpath",
            str(dname_out),
            str(fname_reference),
        ]
    )

if __name__ == "__main__":
    main(
        Path(snakemake.input.fname_reference),
        Path(snakemake.output.fname_insert_bed),
        Path(snakemake.output.dname_out),
        snakemake.params.params,
    )
