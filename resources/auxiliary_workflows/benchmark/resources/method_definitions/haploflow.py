# GROUP: global

"""
Installation from source, no conda enviroment
module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 boost/1.83.0 cmake/3.27.7
. /cluster/work/bewi/members/lfuhrmann/miniconda3/bin/activate
cd /cluster/work/bewi/members/lfuhrmann/
conda activate snakemake
module load eth_proxy
"""

import subprocess
from pathlib import Path
import os

def main(
    fname_bam,
    fname_reference,
    fname_results_snv,
    fname_result_haplos,
    dname_work,
    seq_type,
    threads,
):
    dname_work.mkdir(parents=True, exist_ok=True)

    fname_fastq = str(fname_bam.resolve()).split(".bam")[0] + ".fastq"

    contig_len_filter = 240  # default of Haploflow: 500
    error_rate = 0.0199999996 # default for Illumina
    if seq_type != "illumina":
        error_rate = 0.1 # for long reads

    # execute tool
    subprocess.run(
        [
            "/cluster/work/bewi/members/lfuhrmann/Haploflow/build/haploflow",
            "--read-file",
            fname_fastq,
            "--out",
            dname_work,
            "--log",
            str((dname_work / "log").resolve()),
            "--filter",
            str(contig_len_filter),
            "--error-rate",
            str(error_rate),
            ]
                )

    fname_haplodmf = dname_work / "contigs.fa"
    os.rename(fname_haplodmf.resolve(), fname_result_haplos.resolve())

    # create empty vcf files
    f = open(fname_results_snv, "a")
    f.write("#CHROM POS ID  REF ALT QUAL    FILTER  INFO")
    f.close()


if __name__ == "__main__":
    main(
        Path(snakemake.input.fname_bam),
        Path(snakemake.input.fname_reference),
        Path(snakemake.output.fname_result),
        Path(snakemake.output.fname_result_haplos),
        Path(snakemake.output.dname_work),
        snakemake.wildcards.seq_tech,
        snakemake.threads,
    )
