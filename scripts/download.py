import os
import time
import subprocess
from pathlib import Path

from snakemake.shell import shell


def main(fname_info, accession, restart_times, threads, logfiles):
    outdir = Path(os.path.dirname(fname_info))
    tmpdir = Path(os.path.join(outdir, f"tmp.{accession}"))

    counter = 0
    while True:
        try:
            shell(
                "fasterq-dump --threads {threads} --outdir {outdir} --temp {tmpdir} {accession} > >(tee {logfiles.outfile}) 2>&1"
            )
        except subprocess.CalledProcessError:
            print("Download process crashed, hopefully this is just a fluke...")
            time.sleep(100)

        # make sure the files were actually created
        available_files = list(outdir.glob(f"{accession}*.fastq"))
        if len(available_files) in (1, 2, 3):
            # downloaded SE, PE, varying read number per spot

            # TODO: maybe check all files for >SE reads
            read_len = int(shell.check_output(
                f"bioawk -c fastx '{{{{ bases += length($seq); count++ }}}} END{{{{print int(bases/count)}}}}' {available_files[0]}"
            ).rstrip())

            # TODO: how to get date
            with open(fname_info, "w") as fd:
                fd.write(f"{accession}\t19700101\t{read_len}\n")

            break

        # no files were downloaded, retry...
        shell("echo \"Download failed, restarting\" >> {logfiles.errfile}")
        counter += 1

        if counter > restart_times:
            raise RuntimeError(f"Download failed {counter} times")


if __name__ == "__main__":
    main(
        snakemake.output.fname_info,
        snakemake.wildcards.accession,
        snakemake.params.restart_times,
        snakemake.threads,
        snakemake.log
    )
