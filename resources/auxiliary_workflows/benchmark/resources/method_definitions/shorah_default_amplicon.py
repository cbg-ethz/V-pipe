# GROUP: global
# CONDA: shorah = 1.99.2
# CONDA: biopython = 1.79

import subprocess
from pathlib import Path
from os import listdir
from os.path import isfile, join
from Bio import SeqIO


def main(fname_bam, fname_reference, fname_result, fname_result_haplos, dname_work,seq_tech, genome_size,read_length):
    dname_work.mkdir(parents=True, exist_ok=True)

    if (read_length > genome_size) & (seq_tech == 'illumina'):
        open(fname_result_haplos, 'a').close()
        # create empty vcf files
        f = open(fname_result, 'a')
        f.write("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO")
        f.close()

    else:

        subprocess.run(
            [
                "shorah",
                "amplicon",
                "-b",
                fname_bam.resolve(),
                "-f",
                fname_reference.resolve(),
            ],
            cwd=dname_work,
            check=True,
        )

        (dname_work / "SNVs_0.010000_final.vcf").rename(fname_result)

        mypath = (dname_work).resolve()
        onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
        for file in onlyfiles:
            if file.endswith('reads-support.fas'):
                fname_haplos = (dname_work / "support" / onlyfiles[0]).resolve()
                (dname_work  / file).rename(fname_result_haplos)

        # fix frequency information

        freq_list = []
        for record in SeqIO.parse(fname_result_haplos, "fasta"):
            freq_list.append(float(record.description.split('ave_reads=')[-1]))
        norm_freq_list = [float(i)/sum(freq_list) for i in freq_list]

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
        snakemake.wildcards.seq_tech,
        snakemake.wildcards.genome_size,
        snakemake.wildcards.read_length,
    )
