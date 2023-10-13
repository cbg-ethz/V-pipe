from Bio import SeqIO
from pathlib import Path


def main(fname_reference, fname_insert_bed, genome_size):
    for seq in SeqIO.parse(fname_reference, "fasta"):
        seq_id = seq.id

    with open(fname_insert_bed, "w") as f:
        f.write(
            str(seq_id)
            + "  1       "
            + str(genome_size)
            + "     scheme_INSERT_1 1       +"
        )


if __name__ == "__main__":
    main(
        Path(snakemake.input.fname_reference),
        Path(snakemake.output.fname_insert_bed),
        snakemake.wildcards.genome_size,
    )
