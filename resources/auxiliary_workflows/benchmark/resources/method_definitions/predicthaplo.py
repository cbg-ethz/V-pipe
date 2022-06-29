# GROUP: global
# CONDA: cxx-compiler = 1.4.1
# CONDA: make = 4.3
# CONDA: cmake = 3.22.1
# CONDA: liblapack = 3.9.0
# CONDA: gtest = 1.11.0
# CONDA: samtools = 1.15.1
# CONDA: biopython = 1.79

import subprocess
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def main(fname_bam, fname_reference, fname_result, dname_work):
    dname_work.mkdir(parents=True, exist_ok=True)

    # compile tool
    subprocess.run(
        ["git", "clone", "https://github.com/cbg-ethz/PredictHaplo.git"],
        cwd=dname_work,
        check=True,
    )

    repo_dir = dname_work / "PredictHaplo"
    subprocess.run(
        ["cmake", "-DCMAKE_BUILD_TYPE=Release", "-B", "build", "-S", "."],
        cwd=repo_dir,
        check=True,
    )
    subprocess.run(["cmake", "--build", "build"], cwd=repo_dir, check=True)

    exec_path = repo_dir / "build" / "predicthaplo"

    # prepare environment
    subprocess.run(
        ["samtools", "view", "-h", "-o", dname_work / "reads.sam", fname_bam],
        check=True,
    )

    # execute tool
    ph_prefix = dname_work / "predicthaplo_output"
    subprocess.run(
        [
            exec_path,
            "--sam",
            dname_work / "reads.sam",
            "--reference",
            fname_reference,
            "--prefix",
            ph_prefix,
            "--reconstructed_haplotypes",
            fname_result,
            "--have_true_haplotypes",
            "0",
            "--min_align_score_fraction",
            "-1",
            "--min_overlap_factor",
            "0.1",
            "--entropy_threshold",
            "1e-4",
            "--min_length",
            "1",
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
            Seq(seq),
            id=record.id,
            description=desc,
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
