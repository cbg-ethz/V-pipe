import subprocess
from pathlib import Path

import pandas as pd
import numpy as np
from sklearn import manifold

import seaborn as sns
import matplotlib.pyplot as plt

import editdistance
from Bio import SeqIO


def read_fasta_files(fasta_files):
    tmp = []
    for fname in fasta_files:
        parts = str(fname).split("/")

        if len(parts) == 6:
            _, _, params, _, replicate, _ = parts
            method = None
        elif len(parts) == 7:
            _, _, params, method, _, replicate, _ = parts
        else:
            raise RuntimeError(f"Cannot parse {parts}")

        for record in SeqIO.parse(fname, "fasta"):
            # description actually starts with id
            description = record.description[len(record.id) + 1 :]
            props = dict(pair.split(":") for pair in description.split("|"))

            tmp.append(
                {
                    "method": method,
                    "params": params,
                    "replicate": replicate,
                    "sequence": str(record.seq),
                    "frequency": props.get("freq"),
                }
            )

    return pd.DataFrame(tmp)


def overview_plots(df_haplo, dname_out):
    if df_haplo.empty:
        print("Warning: df_haplo is empty")
        return

    df_haplo["seq_len"] = df_haplo["sequence"].str.len()
    df_long = pd.melt(df_haplo, id_vars=["method", "params", "replicate"]).assign(
        params=lambda x: x["params"].str.replace("_", "\n")
    )
    df_long = df_long[df_long["variable"] != "sequence"]

    g = sns.catplot(
        data=df_long,
        x="params",
        y="value",
        hue="method",
        col="variable",
        kind="box",
        height=10,
    )

    for ax in g.axes.flat:
        ax.tick_params(axis="x", which="major", labelsize=1)

    g.savefig(dname_out / "overview.pdf")


def run_metaquast(predicted_haplos_list, true_haplos_list, workdir):
    df_list = []
    for fname_contigs, fname_truth in zip(predicted_haplos_list, true_haplos_list):
        cwd = workdir / fname_contigs.parent
        # split reference fasta into individual files
        ref_dir = cwd / "haplotype_references"
        ref_dir.mkdir(parents=True, exist_ok=True)

        reference_fname_list = []
        for record in SeqIO.parse(fname_truth, "fasta"):
            fname = ref_dir / f"{record.id}.fasta"
            SeqIO.write(record, fname, "fasta")
            reference_fname_list.append(fname)

        # run quast
        subprocess.run(
            [
                "metaquast",
                "-o",
                cwd,
                "-r",
                ",".join(str(p) for p in reference_fname_list),
                "--min-contig",
                "0",
                "--silent",
                fname_contigs,
            ],
            check=True,
        )

        # parse output
        for res_dir in (cwd / "runs_per_reference").iterdir():
            if res_dir.name.startswith("."):
                continue

            # gather report
            quast_report = pd.read_csv(
                res_dir / "report.tsv",
                sep="\t",
                header=None,
                names=["variable", "value"],
            ).set_index("variable")
            tmp = pd.DataFrame(
                {
                    "contig_count": quast_report.loc["# contigs", "value"],
                    "contig_total_length": quast_report.loc["Total length", "value"],
                    "contig_max_length": quast_report.loc["Largest contig", "value"],
                    "N50": quast_report.loc["N50", "value"],
                    "N75": quast_report.loc["N75", "value"],
                    "L50": quast_report.loc["L50", "value"],
                    "L75": quast_report.loc["L75", "value"],
                },
                index=[0],
            ).astype(int)
            tmp["reference"] = res_dir.name

            # finalize
            _, _, params, method, _, replicate, _ = str(fname_contigs).split("/")
            tmp["params"] = params
            tmp["method"] = method
            tmp["replicate"] = replicate

            df_list.append(tmp)

    return pd.concat(df_list, ignore_index=True)


def plot_quast(df_quast, dname_out):
    for params, group in df_quast.groupby("params"):
        for col in group.select_dtypes(include="number"):
            fig, ax = plt.subplots(figsize=(8, 6))

            sns.boxplot(data=group, y=col, x="method", ax=ax)
            sns.swarmplot(data=group, y=col, x="method", color=".25", ax=ax)

            fig.tight_layout()
            fig.savefig(dname_out / f"{col}__{params}.pdf")


def sequence_embedding(df_pred, df_true, dname_out):
    # compute dissimilarities
    sequence_list = df_pred["sequence"].tolist() + df_true["sequence"].tolist()

    mat = np.zeros(shape=(len(sequence_list), len(sequence_list)))
    print(mat.shape)
    for i, seq1 in enumerate(sequence_list):
        for j, seq2 in enumerate(sequence_list):
            if i >= j:
                continue

            print(i, j)
            mat[i, j] = editdistance.eval(seq1, seq2)

    mat = np.triu(mat) + np.tril(mat.T, 1)  # mirror to make symmetric

    # do MDS
    embedding = manifold.MDS(n_components=2, dissimilarity="precomputed")
    mat_trans = embedding.fit_transform(mat)

    df = pd.concat(
        [
            pd.DataFrame(mat_trans, columns=["MDS0", "MDS1"]),
            pd.concat([df_pred, df_true], axis=0, ignore_index=True),
        ],
        axis=1,
    )
    df["method"] = df["method"].apply(lambda x: "ground_truth" if x is None else x)

    # plot result
    fig, ax = plt.subplots(figsize=(8, 6))

    sns.scatterplot(data=df, x="MDS0", y="MDS1", hue="method", ax=ax)

    fig.savefig(dname_out / "sequence_mds.pdf")


def main(predicted_haplos_list, true_haplos_list, dname_out):
    dname_out.mkdir(parents=True)

    # read data
    df_pred = read_fasta_files(predicted_haplos_list)
    df_true = read_fasta_files(true_haplos_list)

    # create plots
    overview_plots(df_pred, dname_out)

    # quast stuff
    df_quast = run_metaquast(
        predicted_haplos_list, true_haplos_list, dname_out / "quast"
    )
    plot_quast(df_quast, dname_out)

    # MDS
    sequence_embedding(df_pred, df_true, dname_out)


if __name__ == "__main__":
    main(
        [Path(e) for e in snakemake.input.predicted_haplos_list],
        [Path(e) for e in snakemake.input.true_haplos_list],
        Path(snakemake.output.dname_out),
    )
