# GROUP: global
# CONDA: samtools = 1.15.1
# CONDA: biopython = 1.79

"""
To execute PredictHaplo on Euler, I have to load the following:
source /cluster/apps/local/env2lmod.sh # scp doesn't know env2lmod
export MY_MODULEPATH_ROOT=/cluster/work/bewi/nss/modules/
module use $MY_MODULEPATH_ROOT/Core
module load  gcc/5.4.0 openblas/0.2.19
. /cluster/work/bewi/members/lfuhrmann/miniconda3/bin/activate
conda activate snakemake
module load eth_proxy
"""


"""
Note: On Illumina samples, we are receiving this error:
https://github.com/cbg-ethz/PredictHaplo/issues/25

And since there is no newer version to use where single reads are treated
this seems to be not applicalbe for Illumina reads.
"""

import subprocess
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def write_config_file(
    fname_ref, fname_sam, fname_config, ph_prefix, dname_work, fname_result_haplos
):

    read_length = str(fname_sam).split("read_length~")[1].split("__")[0]
    if read_length == "Ten_strain_IAV":
        read_length = 2300
        fname_output_haplos = "predicthaplo_output_global_1_2263.fas"

    with open((dname_work / "config_5V").resolve(), "r") as file:
        content = file.readlines()

    modified_content = content
    modified_content[2] = str(ph_prefix) + "\n"
    modified_content[4] = str(fname_ref) + "\n"
    modified_content[6] = "0" + "\n"
    modified_content[8] = str(fname_sam) + "\n"
    modified_content[10] = "0" + "\n"
    modified_content[12] = str(fname_result_haplos) + "\n"
    modified_content[16] = "10000" + "\n"
    modified_content[20] = "1" + "\n"
    modified_content[22] = str(int(read_length) + 1) + "\n"

    # the following parameters are set such that we don't get the following error:
    # died with <Signals.SIGFPE: 8>.
    modified_content[30] = "-1" + "\n"
    modified_content[34] = "0.1" + "\n"
    modified_content[18] = "0.0001" + "\n"
    modified_content[26] = "1" + "\n"

    with open(fname_config, "w") as file:
        file.writelines(modified_content)


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

    predicthaplo_exe = (
        "/cluster/work/bewi/members/lfuhrmann/PredictHaplo-1.1/PredictHaplo"
    )

    # prepare environment
    subprocess.run(
        ["samtools", "view", "-h", "-o", dname_work / "reads.sam", fname_bam],
        check=True,
    )

    fname_sam = (dname_work / "reads.sam").resolve()

    # copy example config file
    subprocess.run(
        [
            "wget",
            "https://raw.githubusercontent.com/bmda-unibas/PredictHaplo/master/config_5V",
        ],
        cwd=dname_work,
        check=True,
    )

    # modify config file for predicthaplo
    fname_config = (dname_work / "predicthaplo_config").resolve()
    ph_prefix = (dname_work / "predicthaplo_output_").resolve()
    write_config_file(
        Path(fname_reference).resolve(),
        fname_sam,
        fname_config,
        ph_prefix,
        dname_work,
        fname_result_haplos,
    )

    # execute tool
    subprocess.run(
        [
            predicthaplo_exe,
            str(fname_config),
        ],
        check=True,
    )

    # output name t
    read_length = str(fname_sam).split("read_length~")[1].split("__")[0]
    if read_length == "Ten_strain_IAV":
        fname_output_haplos = "predicthaplo_output_global_1_2263.fas"
    else:
        fname_output_haplos = "predicthaplo_output_global_1_5001.fas"

    # clean output file
    record_list = []
    for record in SeqIO.parse((dname_work / fname_output_haplos).resolve(), "fasta"):
        freq = record.description.split("Freq:")[1].split(". Overla")[0]
        record.description = f"freq:{freq}"
        seq_nodel = str(record.seq).replace("-", "")
        record.seq = Seq(seq_nodel)
        record_list.append(record)

    SeqIO.write(record_list, fname_result_haplos, "fasta")

    # create empty vcf files
    f = open(fname_results_snv, "a")
    f.write("#CHROM     POS     ID      REF     ALT     QUAL    FILTER  INFO")
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
