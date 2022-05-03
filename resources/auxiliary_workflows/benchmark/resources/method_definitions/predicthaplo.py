# GROUP: global
# CONDA: cxx-compiler = 1.4.1
# CONDA: cmake = 3.22.1
# CONDA: liblapack = 3.9.0
# CONDA: gtest = 1.11.0
# CONDA: samtools = 1.15.1
# CONDA: biopython = 1.79

import subprocess
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def main(fname_bam, fname_reference, fname_result, dname_work):
    dname_work.mkdir(parents=True, exist_ok=True)

    # compile tool
    subprocess.run(
        ["git", "clone", "https://github.com/cbg-ethz/PredictHaplo"],
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
            "--have_true_haplotypes",
            "0",
            "--min_align_score_fraction",
            "-1",
            "--min_overlap_factor",
            "0.1",
        ],
        check=True,
    )

    # aggregate result
    record_list = []
    for path in dname_work.glob(f"{ph_prefix.name}global*.fas"):
        for record in SeqIO.parse(path, "fasta"):
            seq = record.seq.split("EndOfComments")[1]
            window = path.name[len(f"{ph_prefix.name}global_") : -len(".fas")]
            id_ = f"{window}__{record.id}"

            rec = SeqRecord(
                seq,
                id=id_,
                name=id_,
                description=str(record.seq).split(";")[1],
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
