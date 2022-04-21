# GROUP: global
# CONDA: haploclique = 1.3.1

import subprocess
from pathlib import Path


def main(fname_bam, fname_reference, fname_marker, dname_work):
    dname_work.mkdir(parents=True, exist_ok=True)

    subprocess.run(
        [
            "haploclique",
            "--no_singletons",
            "--no_prob0",
            "--edge_quasi_cutoff_cliques=0.85",
            "--edge_quasi_cutoff_mixed=0.85",
            "--edge_quasi_cutoff_single=0.8",
            "--min_overlap_cliques=0.6",
            "--min_overlap_single=0.5",
            "--limit_clique_size=3",
            "--max_cliques=10000",
            fname_bam.resolve(),
        ],
        cwd=dname_work,
        check=True,
    )

    (dname_work / "quasispecies.fasta").rename(fname_marker.parent / "haplotypes.fasta")


if __name__ == "__main__":
    main(
        Path(snakemake.input.fname_bam),
        Path(snakemake.input.fname_reference),
        Path(snakemake.output.fname_marker),
        Path(snakemake.output.dname_work),
    )
