# GROUP: global
# CONDA: cxx-compiler = 1.4.1
# CONDA: cmake = 3.22.1
# CONDA: liblapack = 3.9.0
# CONDA: gtest = 1.11.0

import subprocess
from pathlib import Path


def main(fname_bam, fname_reference, fname_marker, dname_work):
    dname_work.mkdir(parents=True, exist_ok=True)

    # compile tool
    subprocess.run(
        ["git", "clone", "https://github.com/cbg-ethz/PredictHaplo"], cwd=dname_work
    )

    repo_dir = dname_work / "PredictHaplo"
    subprocess.run(
        ["cmake", "-DCMAKE_BUILD_TYPE=Release", "-B", "build", "-S", "."], cwd=repo_dir
    )
    subprocess.run(["cmake", "--build", "build"], cwd=repo_dir)

    exec_path = repo_dir / "build" / "predicthaplo"

    # prepare environment
    subprocess.run(
        ["samtools", "view", "-h", "-o", dname_work / "reads.sam", fname_bam]
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
        ]
    )

    # TODO: choose correct results file
    (dname_work / "TODO").rename(fname_marker.parent / "haplotypes.fasta")


if __name__ == "__main__":
    main(
        Path(snakemake.input.fname_bam),
        Path(snakemake.input.fname_reference),
        Path(snakemake.output.fname_marker),
        Path(snakemake.output.dname_work),
    )
