# GROUP: local
# CONDA: cliquesnv = 2.0.3
# CONDA: samtools = 1.15.1
# CONDA: biopython = 1.79

import subprocess
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


def main(fname_bam, fname_reference, fname_results_snv, fname_result_haplos, dname_work, seq_type, threads=1):
    dname_work.mkdir(parents=True, exist_ok=True)

    # prepare environment
    subprocess.run(
        ["samtools", "view", "-h", "-o", dname_work / "reads.sam", fname_bam],
        check=True,
    )

    cliquesnv_mode = None
    if seq_type == "illumina":
        cliquesnv_mode = "snv-illumina-vc"
    elif seq_type == "pacbio":
        cliquesnv_mode = "snv-pacbio-vc"
    elif seq_type == "nanopore":
        cliquesnv_mode = "snv-pacbio-vc"
    else:
        raise RuntimeError(f"Invalid sequence technology: {seq_type}")

    # execute tool
    subprocess.run(
        [
            "cliquesnv",
            "-m",
            cliquesnv_mode,
            "-in",
            dname_work / "reads.sam",
            "-outDir",
            dname_work / "output",
            "-tf", # parameter to detect low-frequent mutations
            "0.01",
            "-tl", # maximal runtime parameter
            "2*428400", # 2*119h*60*60
            "-threads",
            str(threads),
            "-Xms612m",
            "-Xmx40g",
        ],
        check=True,
    )
    (dname_work / "output" / "reads.vcf").rename(fname_results_snv)

    # create empty haplotype files
    open(fname_result_haplos, 'a').close()


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
