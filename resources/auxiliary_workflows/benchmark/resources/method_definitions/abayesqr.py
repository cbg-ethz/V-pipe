# GROUP: global
# CONDA: abayesqr = 1.0
# CONDA: samtools = 1.15.1

import subprocess
from pathlib import Path


def main(fname_bam, fname_reference, fname_marker, dname_work):
    dname_work.mkdir(parents=True, exist_ok=True)

    # create config
    fname_sam = (dname_work / "reads.sam").resolve()
    config = f"""filename of reference sequence (FASTA) : {fname_reference.resolve()}
filname of the aligned reads (sam format) : {fname_sam}
paired-end (1 = true, 0 = false) : 1
SNV_thres : 0.01
reconstruction_start : 1
reconstruction_stop: 1300
min_mapping_qual : 60
min_read_length : 150
max_insert_length : 250
characteristic zone name : test
seq_err (assumed sequencing error rate(%)) : 0.1
MEC improvement threshold : 0.0395"""
    (dname_work / "config.txt").write_text(config)

    # execute tool
    subprocess.run(
        ["samtools", "view", "-h", "-o", fname_sam, fname_bam],
        check=True,
    )

    subprocess.run(["aBayesQR", "config.txt"], cwd=dname_work, check=True)

    # aggregate results
    # TODO (<zone name>_ViralSeq.txt)
    # (dname_work / "output" / "quasispecies.fasta").rename(
    #     fname_marker.parent / "haplotypes.fasta"
    # )


if __name__ == "__main__":
    main(
        Path(snakemake.input.fname_bam),
        Path(snakemake.input.fname_reference),
        Path(snakemake.output.fname_marker),
        Path(snakemake.output.dname_work),
    )
