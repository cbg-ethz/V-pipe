# GROUP: global
# CONDA: quasirecomb = 1.2
# CONDA: biopython = 1.79

import subprocess
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def main(fname_bam, fname_reference, fname_result, dname_work):
    dname_work.mkdir(parents=True, exist_ok=True)

    subprocess.run(
        ["quasirecomb", "-i", fname_bam, "-o", dname_work / "output"], check=True
    )
    fname_quasi = dname_work / "output" / "quasispecies.fasta"

    # clean output file
    record_list = []
    for record in SeqIO.parse(fname_quasi, "fasta"):
        id_, freq = record.id.split("_")

        rec = SeqRecord(
            record.seq,
            id=id_,
            description=f"freq:{float(freq)}",
        )
        record_list.append(rec)

    SeqIO.write(record_list, fname_result, "fasta")


if __name__ == "__main__":
    main(
        Path(snakemake.input.fname_bam),
        Path(snakemake.input.fname_reference),
        Path(snakemake.output.fname_result),
        Path(snakemake.output.dname_work),
    )
