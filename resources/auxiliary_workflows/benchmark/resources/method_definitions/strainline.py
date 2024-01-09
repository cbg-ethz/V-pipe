# GROUP: global
# CONDA: minimap2 = 2.24
# CONDA: spoa
# CONDA: samtools
# CONDA: dazz_db
# CONDA: daligner
# CONDA: metabat2
# CONDA: biopython = 1.79

import subprocess
from pathlib import Path
from Bio import SeqIO


def fastq_to_fasta(input_file, output_file):
    """
    Converts a Fastq file to a Fasta file.

    Args:
        input_file (str): Path to the input Fastq file.
        output_file (str): Path to save the output Fasta file.
    """

    with open(output_file, "w") as f_out:
        for record in SeqIO.parse(input_file, "fastq"):
            SeqIO.write(record, f_out, "fasta")



def main(fname_bam,
         fname_reference,
         fname_result,
         fname_result_haplos,
         dname_work,
         seq_type,
         threads,):

    if seq_type == "illumina":
        # strainline doesn't work for illumina
        # create fake files such that snakemake is happy

        # create empty haplotype files
        open(fname_result_haplos, "a").close()

        # create empty vcf files
        f = open(fname_results_snv, "a")
        f.write("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO")
        f.close()
    else:
        # run strainline

        if seq_type == "pacbio":
            seq_platform = "pb"
        elif seq_type == "nanopore":
            seq_platform = "ont"

        # fastq to fasta
        fname_reads_fastq = str(fname_bam.resolve()).split("bam")[0]+"fastq"
        fname_reads_fasta = str(fname_bam.resolve()).split("bam")[0]+"fasta"
        fastq_to_fasta(fname_reads_fastq, fname_reads_fasta)

        dname_work.mkdir(parents=True, exist_ok=True)

        # execute tool
        subprocess.run(
            [
                "strainline.sh",
                "-i",
                fname_reads_fasta,
                "-o",
                str(dname_work),
                "-p",
                seq_platform,
                "-threads",
                str(threads),
            ],
            check=True,
        )





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
