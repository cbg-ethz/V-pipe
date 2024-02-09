# GROUP: global
# CONDA: cxx-compiler = 1.4.1
# CONDA: make = 4.3
# CONDA: cmake = 3.22.1
# CONDA: rhash = 1.4.4
# CONDA: liblapack = 3.9.0
# CONDA: gtest = 1.11.0
# CONDA: samtools = 1.15.1
# CONDA: biopython = 1.79

"""
To execute PredictHaplo on Euler, I have to load the following:
source /cluster/apps/local/env2lmod.sh # scp doesn't know env2lmod
export MY_MODULEPATH_ROOT=/cluster/work/bewi/nss/modules/
module use $MY_MODULEPATH_ROOT/Core
module load  gcc/5.4.0 openblas/0.2.19
. /cluster/work/bewi/members/lfuhrmann/miniconda3/bin/activate
cd /cluster/work/bewi/members/lfuhrmann/
conda activate snakemake
"""


import subprocess
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def write_config_file(fname_ref, fname_sam,fname_config, ph_prefix):
    with open('original_file.txt', 'r') as file:
        content = file.readlines()

    modified_content =  content
    modified_content[2] = str(ph_prefix)
    modified_content[4] = str(fname_ref)
    modified_content[6] = "0"
    modified_content[8] = str(fname_sam)
    modified_content[10] = "0"
    modified_content[12] = str(fname_ref)
    modified_content[16] = "10000"

    with open(fname_config, 'w') as file:
        file.writelines(modified_content)


def main(fname_bam,
         fname_reference,
         fname_results_snv,
         fname_result_haplos,
         dname_work,
         seq_type,
         threads,
):
    dname_work.mkdir(parents=True, exist_ok=True)

    predicthaplo_exe = "/cluster/work/bewi/members/lfuhrmann/PredictHaplo-1.1/PredictHaplo"

    # prepare environment
    subprocess.run(
        ["samtools", "view", "-h", "-o", dname_work / "reads.sam", fname_bam],
        check=True,
    )

    fname_sam = (dname_work / "reads.sam").resolve()

    # copy example config file
    subprocess.run(
        ["wget", "https://github.com/bmda-unibas/PredictHaplo/blob/master/config_5V"],
        cwd=dname_work,
        check=True,
    )

    # modify config file for predicthaplo
    fname_config = (dname_work / "predicthaplo_config").resolve()
    ph_prefix = (dname_work / "predicthaplo_output").resolve()
    write_config_file(Path(fname_reference).resolve(), fname_sam, fname_config, ph_prefix)

    # execute tool
    subprocess.run(
        [
            predicthaplo_exe,
            str(fname_config),
        ],
        check=True,
    )

    # clean output file
    record_list = []
    for record in SeqIO.parse(fname_result, "fasta"):
        props = str(record.seq).split("EndOfComments")[0]
        seq = str(record.seq).split("EndOfComments")[-1]

        desc = props.split(";")[1]  # only keep frequency

        rec = SeqRecord(
            Seq(seq.replace("-", "")),
            id=record.id,
            description=desc,
        )
        record_list.append(rec)

    SeqIO.write(record_list, fname_result, "fasta")

    # create empty vcf files
    f = open(fname_results_snv, "a")
    f.write("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO")
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
