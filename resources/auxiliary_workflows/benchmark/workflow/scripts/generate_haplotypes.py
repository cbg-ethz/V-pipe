from pathlib import Path

import numpy as np
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import seaborn as sns
import matplotlib.pylab as plt


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


def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet.
        N: number of colors.

    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    """

    if type(cmap) == str:
        cmap = get_cmap(cmap)
    colors_i = concatenate((linspace(0, 0.2, N), (0.0, 0.0, 0.0, 0.0)))
    colors_rgba = cmap(colors_i)
    indices = linspace(0, 0.2, N + 1)
    cdict = {}
    for ki, key in enumerate(("red", "green", "blue")):
        cdict[key] = [
            (indices[i], colors_rgba[i - 1, ki], colors_rgba[i, ki])
            for i in xrange(N + 1)
        ]
    # Return colormap object.
    return matplotlib.colors.LinearSegmentedColormap(cmap.name + "_%d" % N, cdict, 1024)


def hamming(chaine1, chaine2):
    return sum(c1 != c2 for c1, c2 in zip(chaine1, chaine2))


def compute_pairwise_hamming(haplos_list, haplo_ids):

    n_haplos = len(haplos_list)
    dist = np.zeros((n_haplos, n_haplos))
    for i_seq, seq in enumerate(haplos_list):
        for j_seq, seq2 in enumerate(haplos_list):
            dist[i_seq][j_seq] = hamming(seq, seq2)

    distance_matrix1 = pd.DataFrame(columns=haplo_ids, index=haplo_ids, data=dist)
    return distance_matrix1


def mutate(master_seq, distance):
    seq_haplotype = np.asarray(list(master_seq))
    ground_truth = {"type": [], "position": [], "variant": []}

    position_list = np.random.choice(
        np.arange(len(master_seq)), size=distance, replace=False
    )

    base_set = set(BASE_LIST)
    for pos in position_list:
        # make sure we mutate to base different from reference
        cur_base_list = list(base_set - {seq_haplotype[pos]})
        seq_haplotype[pos] = np.random.choice(cur_base_list)

    ground_truth["type"].extend(["mutation"] * position_list.shape[0])
    ground_truth["position"].extend(position_list)
    ground_truth["variant"].extend(seq_haplotype[position_list])

    return "".join(seq_haplotype), pd.DataFrame(ground_truth)


def generate_haplotype_groups(
    master_seq, n_group1, n_group2, d_group1, d_group2, d_group12
):
    """d_group1, d_group2, d_group12: must be divisible by 2 (implementation reasons)."""
    genome_length = len(master_seq)

    # generate master of group 1 and group 2
    master_group1, _ = mutate(master_seq, int(d_group12 / 2))
    master_group2, _ = mutate(master_seq, int(d_group12 / 2))

    ground_truth_list = []
    haplo_list = []

    for i in range(n_group1):
        seq, ground_truth = mutate(master_group1, int(d_group1 / 2))
        haplo_list.append(seq)
        ground_truth_list.append(ground_truth)

    for i in range(n_group2):
        seq, ground_truth = mutate(master_group2, int(d_group2 / 2))
        haplo_list.append(seq)
        ground_truth_list.append(ground_truth)

    return haplo_list, ground_truth_list


def plot_pairwise_distances(distances, dname_work):
    mask = np.zeros_like(distances)
    mask[np.triu_indices_from(mask)] = True
    cmap = sns.cubehelix_palette(start=2.8, rot=0.1, light=0.9, n_colors=8)
    plt.figure(figsize=(45, 10))
    with sns.axes_style("white"):
        plt_sns = sns.heatmap(
            distances, mask=mask, square=True, annot=True, cbar=False, cmap=cmap
        )
        figure = plt_sns.get_figure()
        figure.suptitle("Pairwise hamming distances", fontsize=16)
        figure.savefig(str(dname_work) + "/local_haplo_pairwise_distance.pdf", dpi=400)


def main(
    fname_reference,
    fname_groundtruth,
    fname_fasta,
    dname_work,
    haplotype_generation,
    params,
):
    """Create master sequence, infer haplotypes and simulate reads."""
    # initial setup
    # np.random.seed(42)
    dname_work.mkdir(parents=True, exist_ok=True)

    # generate random master sequence
    master_name = "MasterSequence"
    seq_master = "".join(np.random.choice(BASE_LIST, size=params["genome_size"]))
    fname_reference.write_text(f">{master_name}\n{seq_master}\n")

    if haplotype_generation == "distance":
        (
            n_group1,
            n_group2,
            d_group12,
            d_group1,
            d_group2,
            freq_dist,
            freq_param,
        ) = params["haplos"].split("@")
        n_group1 = int(n_group1)
        n_group2 = int(n_group2)
        d_group12 = int(d_group12)
        d_group1 = int(d_group1)
        d_group2 = int(d_group2)
        n_haplo = int(n_group1) + int(n_group2)
        # obtain haplotype frequencies
        if freq_dist == "dirichlet":
            if type(freq_param) == str:
                alpha = [float(freq) for freq in freq_param.split(":")]
            if type(freq_param) == np.float64:
                if np.isnan(freq_param):
                    alpha = np.ones(n_haplo)

            freq_list = np.random.dirichlet(alpha)

        elif freq_dist == "geom":
            freq_param = float(freq_param)
            freq_list = np.asarray([freq_param ** (i + 1) for i in range(n_haplo)])
            freq_list = freq_list / np.sum(freq_list)

        ground_truth_list = []
        filelist_sam = []
        filelist_fastq = []

        haplo_list, ground_truth_list_temp = generate_haplotype_groups(
            seq_master,
            n_group1,
            n_group2,
            d_group1,
            d_group2,
            d_group12,
        )

        haplo_ids = [f"haplotype{i:04}" for i in range(n_haplo)]

        # compute pairwise distances
        distances = compute_pairwise_hamming(haplo_list, haplo_ids)
        distances.to_csv(str(dname_work) + "/pairwise_hamming_distance.csv")

        # plot heatmap of pairwise distances
        plot_pairwise_distances(distances, dname_work)

    elif haplotype_generation == "mutation_rate":
        mutation_rate, insertion_rate, deletion_rate, haplotype_pattern = params[
            "haplos"
        ].split("@")
        mutation_rate = float(mutation_rate)
        insertion_rate = float(insertion_rate)
        deletion_rate = float(deletion_rate)
        # infer haplotype sequences
        freq_list = [float(freq) for freq in haplotype_pattern.split(":")]
        assert (
            sum(freq_list) == 1
        ), f"Invalid haplotype pattern: {haplotype_pattern}, sum is {sum(freq_list)}"

        n_haplo = len(freq_list)

        ground_truth_list = []
        ground_truth_list_temp = []
        haplo_list = []
        for i, freq in enumerate(freq_list):
            # generate haplotype
            haplotype_name = f"haplotype{i:04}"
            seq_haplotype, ground_truth = generate_haplotype(
                seq_master,
                mutation_rate,
                insertion_rate,
                deletion_rate,
            )
            ground_truth_list_temp.append(ground_truth)
            haplo_list.append(seq_haplotype)

        haplo_ids = [f"haplotype{i:04}" for i in range(n_haplo)]

    record_list = []
    for i, freq in enumerate(freq_list):
        # generate haplotype
        ground_truth = ground_truth_list_temp[i]
        seq_haplotype = haplo_list[i]
        haplotype_name = haplo_ids[i]

        ground_truth["haplotype"] = haplotype_name

        # save haplotype in FASTA
        fname_haplotype = dname_work / f"{haplotype_name}.fasta"
        fname_haplotype.write_text(f">{haplotype_name}\n{seq_haplotype}\n")

        rec = SeqRecord(
            Seq(seq_haplotype),
            id=haplotype_name,
            description=f"freq:{freq}",
        )
        record_list.append(rec)

        ground_truth_list.append(ground_truth)

    pd.concat(ground_truth_list).to_csv(fname_groundtruth)
    SeqIO.write(record_list, fname_fasta, "fasta")


if __name__ == "__main__":
    main(
        Path(snakemake.output.fname_reference),
        Path(snakemake.output.fname_groundtruth),
        Path(snakemake.output.fname_fasta),
        Path(snakemake.output.dname_work),
        snakemake.params.haplotype_generation,
        snakemake.params.params,
    )
