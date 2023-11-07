# GROUP: global
# CONDA: libshorah
# CONDA: biopython = 1.79
# PIP: git+https://github.com/LaraFuhrmann/shorah@master

import subprocess
from pathlib import Path
from os import listdir
from os.path import isfile, join
from Bio import SeqIO
import gzip


def gunzip(source_filepath, dest_filepath, block_size=65536):
    with gzip.open(source_filepath, "rb") as s_file, open(
        dest_filepath, "wb"
    ) as d_file:
        while True:
            block = s_file.read(block_size)
            if not block:
                break
            else:
                d_file.write(block)


def main(
    fname_bam,
    fname_reference,
    fname_insert_bed,
    fname_results_snv,
    fname_result_haplos,
    dname_work,
    threads,
):
    genome_size = str(fname_bam).split("genome_size~")[1].split("__coverage")[0]
    alpha = 0.000001
    n_max_haplotypes = 500
    n_mfa_starts = 1
    win_min_ext = 0.85

    read_length = str(fname_bam).split("read_length~")[1].split("__")[0]
    if read_length == "Ten_strain_IAV":
        sampler = "learn_error_params"
        win_min_ext = 0.5
    else:
        sampler = "use_quality_scores"

    dname_work.mkdir(parents=True, exist_ok=True)
    if fname_insert_bed == []:
        subprocess.run(
            [
                "shorah",
                "shotgun",
                "-b",
                fname_bam.resolve(),
                "-f",
                Path(fname_reference).resolve(),
                "--mode",
                str(sampler),
                "--alpha",
                str(alpha),
                "--n_max_haplotypes",
                str(n_max_haplotypes),
                "--n_mfa_starts",
                str(n_mfa_starts),
                "--win_min_ext",
                str(win_min_ext),
                "--threads",
                str(threads),
                "--extended_window_mode",
            ],
            cwd=dname_work,
        )
    else:
        # insert bed file is there
        subprocess.run(
            [
                "shorah",
                "shotgun",
                "-b",
                fname_bam.resolve(),
                "-f",
                Path(fname_reference).resolve(),
                "-z",
                Path(fname_insert_bed).resolve(),
                "--mode",
                str(sampler),
                "--alpha",
                str(alpha),
                "--n_max_haplotypes",
                str(n_max_haplotypes),
                "--n_mfa_starts",
                str(n_mfa_starts),
                "--win_min_ext",
                str(win_min_ext),
                "--threads",
                str(threads),
                "--extended_window_mode",
            ],
            cwd=dname_work,
        )

    (dname_work / "snv" / "SNVs_0.010000_final.vcf").rename(fname_results_snv)
    open(fname_result_haplos, "a").close()


if __name__ == "__main__":
    main(
        Path(snakemake.input.fname_bam),
        Path(snakemake.input.fname_reference),
        snakemake.input.fname_insert_bed,
        Path(snakemake.output.fname_result),
        Path(snakemake.output.fname_result_haplos),
        Path(snakemake.output.dname_work),
        snakemake.threads,
    )
