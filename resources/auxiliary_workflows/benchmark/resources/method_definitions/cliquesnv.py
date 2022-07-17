# GROUP: global
# CONDA: cliquesnv = 2.0.3
# CONDA: samtools = 1.15.1
# CONDA: biopython = 1.79

import subprocess
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


def main(fname_bam, fname_reference, fname_result, dname_work, seq_type, thread_num):
    dname_work.mkdir(parents=True, exist_ok=True)

    # prepare environment
    subprocess.run(
        ["samtools", "view", "-h", "-o", dname_work / "reads.sam", fname_bam],
        check=True,
    )

    cliquesnv_mode = None
    if seq_type == "illumina":
        cliquesnv_mode = "snv-illumina"
    elif seq_type == "pacbio":
        cliquesnv_mode = "snv-pacbio"
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
            "-threads",
            str(thread_num),
            "-Xms512m",
            "-Xmx20g",
        ],
        check=True,
    )
    fname_cliquesnv = dname_work / "output" / "reads.fasta"

    # fix frequency information
    if fname_cliquesnv.exists():
        record_list = []
        for record in SeqIO.parse(fname_cliquesnv, "fasta"):
            _, _, freq = record.id.split("_")
            record.description = f"freq:{freq}"

            seq_nodel = str(record.seq).replace("-", "")
            record.seq = Seq(seq_nodel)

            record_list.append(record)
    else:
        record_list = []

    SeqIO.write(record_list, fname_result, "fasta")


if __name__ == "__main__":
    main(
        Path(snakemake.input.fname_bam),
        Path(snakemake.input.fname_reference),
        Path(snakemake.output.fname_result),
        Path(snakemake.output.dname_work),
        snakemake.wildcards.seq_tech,
        snakemake.threads,
    )
