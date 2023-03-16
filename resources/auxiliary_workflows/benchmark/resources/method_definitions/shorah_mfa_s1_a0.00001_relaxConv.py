# GROUP: global
# CONDA: boost = 1.77.0
# CONDA: htslib = 1.14
# CONDA: biopython = 1.79
# PIP: git+https://github.com/LaraFuhrmann/shorah@feature-new-inference


"""
Configuration:
alpha = 0.00001
relaxConv
start 1 (not so random initialization, so one start should be sufficient)

--> purpose: compare if smaller alpha is better.

"""

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
    fname_bam, fname_reference, fname_results_snv, fname_result_haplos, dname_work
):
    genome_size = str(fname_bam).split("genome_size~")[1].split("__coverage")[0]
    alpha = 0.00001
    n_max_haplotypes = 100
    n_mfa_starts = 1

    dname_work.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "shorah",
            "shotgun",
            "-b",
            fname_bam.resolve(),
            "-f",
            fname_reference.resolve(),
            "-w",
            str(genome_size),
            "--inference",
            "mean_field_approximation",
            "--alpha",
            str(alpha),
            "--n_max_haplotypes",
            str(n_max_haplotypes),
            "--n_mfa_starts",
            str(n_mfa_starts),
        ],
        cwd=dname_work,
    )

    (dname_work / "snv" / "SNVs_0.010000_final.vcf").rename(fname_results_snv)

    mypath = (dname_work / "support").resolve()
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    print("onlyfiles", onlyfiles)
    for file in onlyfiles:
        if "reads-support.fas" in file:
            file_name = onlyfiles[0]
            fname_haplos = (dname_work / "support" / onlyfiles[0]).resolve()
            if file.endswith(".gz"):
                fname_zipped = (dname_work / "support" / onlyfiles[0]).resolve()
                fname_haplos = onlyfiles[0].split(".gz")[0]
                fname_unzipped = (dname_work / "support" / fname_haplos).resolve()
                # unzip
                gunzip(fname_zipped, fname_result_haplos)

            elif file.endswith(".fas"):
                fname_haplos = (dname_work / "support" / onlyfiles[0]).resolve()
                (dname_work / "support" / file).rename(fname_result_haplos)

    # fix frequency information

    freq_list = []
    for record in SeqIO.parse(fname_result_haplos, "fasta"):
        freq_list.append(float(record.description.split("ave_reads=")[-1]))
    norm_freq_list = [float(i) / sum(freq_list) for i in freq_list]

    record_list = []
    for idx, record in enumerate(SeqIO.parse(fname_result_haplos, "fasta")):
        record.description = f"freq:{norm_freq_list[idx]}"
        record_list.append(record)
    SeqIO.write(record_list, fname_result_haplos, "fasta")


if __name__ == "__main__":
    main(
        Path(snakemake.input.fname_bam),
        Path(snakemake.input.fname_reference),
        Path(snakemake.output.fname_result),
        Path(snakemake.output.fname_result_haplos),
        Path(snakemake.output.dname_work),
    )
