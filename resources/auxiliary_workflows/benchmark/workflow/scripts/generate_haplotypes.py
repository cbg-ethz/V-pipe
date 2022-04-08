import shutil
import tempfile
import fileinput
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd


BASE_LIST = list("TCGA")


def generate_haplotype(seq_master, mutation_rate=0, insertion_rate=0, deletion_rate=0):
    """Generate haplotype from master sequence."""
    # TODO: allow indels of length >1

    seq_haplotype = np.asarray(list(seq_master))
    ground_truth = {"type": [], "position": [], "variant": []}

    # deletions
    deletion_count = int(len(seq_haplotype) * deletion_rate)
    position_list = np.random.choice(
        np.arange(len(seq_haplotype)), size=deletion_count, replace=False
    )
    seq_haplotype = np.delete(seq_haplotype, position_list)

    ground_truth["type"].extend(["deletion"] * position_list.shape[0])
    ground_truth["position"].extend(position_list)
    ground_truth["variant"].extend(["-"] * position_list.shape[0])

    # mutations
    mutation_count = int(len(seq_haplotype) * mutation_rate)
    position_list = np.random.choice(
        np.arange(len(seq_haplotype)), size=mutation_count, replace=False
    )

    base_set = set(BASE_LIST)
    for pos in position_list:
        # make sure we mutate to base different from reference
        cur_base_list = list(base_set - {seq_haplotype[pos]})
        seq_haplotype[pos] = np.random.choice(cur_base_list)

    ground_truth["type"].extend(["mutation"] * position_list.shape[0])
    ground_truth["position"].extend(position_list)
    ground_truth["variant"].extend(seq_haplotype[position_list])

    # insertions
    insertion_count = int(len(seq_haplotype) * insertion_rate)
    position_list = np.random.choice(
        np.arange(len(seq_haplotype)), size=insertion_count, replace=False
    )

    insertion_list = np.random.choice(BASE_LIST, size=insertion_count)
    seq_haplotype = np.insert(
        seq_haplotype,
        position_list,
        insertion_list,
    )

    ground_truth["type"].extend(["insertion"] * position_list.shape[0])
    ground_truth["position"].extend(position_list)
    ground_truth["variant"].extend(insertion_list)

    return "".join(seq_haplotype), pd.DataFrame(ground_truth)



def main(
        fname_reference, fname_groundtruth, dname_work, params
):
    """Create master sequence, infer haplotypes and simulate reads."""
    # initial setup
    # np.random.seed(42)
    dname_work.mkdir(parents=True, exist_ok=True)

    # generate random master sequence
    master_name = "MasterSequence"
    seq_master = "".join(np.random.choice(BASE_LIST, size=params["genome_size"]))
    fname_reference.write_text(f">{master_name}\n{seq_master}\n")

    # infer haplotype sequences
    freq_list = [float(freq) for freq in params["haplotype_pattern"].split(":")]
    assert (
        sum(freq_list) == 1
    ), f"Invalid haplotype pattern: {params['haplotype_pattern']}"

    ground_truth_list = []
    for i, freq in enumerate(freq_list):
        # generate haplotype
        haplotype_name = f"haplotype{i:04}"
        seq_haplotype, ground_truth = generate_haplotype(
            seq_master,
            params["mutation_rate"],
            params["insertion_rate"],
            params["deletion_rate"],
        )

        ground_truth["haplotype"] = haplotype_name

        coverage_haplotype = int(params["coverage"] * freq)

        # save haplotype in FASTA
        fname_haplotype = dname_work / f"{haplotype_name}.fasta"
        fname_haplotype.write_text(f">{haplotype_name}\n{seq_haplotype}\n")

        ground_truth_list.append(ground_truth)

    pd.concat(ground_truth_list).to_csv(fname_groundtruth)


if __name__ == "__main__":
    main(
        Path(snakemake.output.fname_reference),
        Path(snakemake.output.fname_groundtruth),
        Path(snakemake.output.dname_work),
        snakemake.params.params,
    )
