# GROUP: global
# CONDA: libshorah
# CONDA: fuc
# CONDA: biopython = 1.79
# PIP: pandas
# PIP: git+https://github.com/cbg-ethz/VILOCA@master

import subprocess
from pathlib import Path
from os import listdir
from os.path import isfile, join
from Bio import SeqIO
import gzip
from fuc import pyvcf
import pandas as pd
import shutil


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


def filter_snvs(fname_vcf_viloca, fname_results_snv, posterior_threshold):
    vf = pyvcf.VcfFrame.from_file(fname_vcf_viloca)
    info_strings = (
        '{"'
        + vf.df.INFO.str.split(";")
        .str.join('","')
        .str.replace("=", '":"')
        .str.replace('"",', "")
        + '"}'
    )
    info_df = pd.json_normalize(info_strings.apply(eval))
    filter_out = info_df["Post1"].astype(float) > posterior_threshold
    vf.df = vf.df[filter_out]
    vf.to_file(fname_results_snv)


def main(
    fname_bam,
    fname_reference,
    fname_insert_bed,
    fname_results_snv,
    fname_result_haplos,
    dname_work,
    threads=1,
):
    genome_size = str(fname_bam).split("genome_size~")[1].split("__coverage")[0]
    alpha = 0.00001
    n_max_haplotypes = 500
    n_mfa_starts = 1
    win_min_ext = 0.85
    viloca_posterior_threshold = 0.0

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
                "viloca",
                "run",
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
                "--threshold",
                str(viloca_posterior_threshold),
            ],
            cwd=dname_work,
        )
    else:
        # insert bed file is there
        subprocess.run(
            [
                "viloca",
                "run",
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
                "--threshold",
                str(viloca_posterior_threshold),
            ],
            cwd=dname_work,
        )

    # here are all snvs regardless their posterior: (dname_work / "snv" / "SNVs_0.010000_final.vcf")

    # filter out snvs with low posterior, threshold = 0.9 as default in shorah
    fname_vcf_viloca = str(
        Path(dname_work / "snv" / "SNVs_0.010000_final.vcf").resolve()
    )
    filter_snvs(
        fname_vcf_viloca,
        str(Path(fname_results_snv).resolve()),
        posterior_threshold=0.9,
    )

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
                shutil.copy(
                    (dname_work / "support" / file).resolve(), fname_result_haplos
                )

    # fix frequency information

    freq_list = []
    post_list = []
    for record in SeqIO.parse(fname_result_haplos, "fasta"):
        freq_list.append(float(record.description.split("ave_reads=")[-1]))
        post_list.append(
            float(record.description.split("posterior=")[-1].split(" ave_")[0])
        )
    norm_freq_list = [float(i) / sum(freq_list) for i in freq_list]

    record_list = []
    for idx, record in enumerate(SeqIO.parse(fname_result_haplos, "fasta")):
        record.description = f"posterior:{post_list[idx]}|freq:{norm_freq_list[idx]}"
        record_list.append(record)
    SeqIO.write(record_list, fname_result_haplos, "fasta")


if __name__ == "__main__":
    main(
        Path(snakemake.input.fname_bam),
        Path(snakemake.input.fname_reference),
        snakemake.input.fname_insert_bed,
        Path(snakemake.output.fname_result),
        Path(snakemake.output.fname_result_haplos),
        Path(snakemake.output.dname_work),
    )
