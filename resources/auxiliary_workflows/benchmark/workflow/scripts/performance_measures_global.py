import functools
import subprocess
from pathlib import Path

import pandas as pd
import numpy as np
from sklearn import manifold

import seaborn as sns
import matplotlib.pyplot as plt

import editdistance
from Bio import SeqIO

from tqdm import tqdm


def read_fasta_files(fasta_files, with_method=True):
    tmp = []
    for fname in tqdm(fasta_files, desc="Read FASTA files"):
        parts = str(fname).split("/")

        if with_method:
            params = parts[-5]
            method = parts[-4]
        else:
            params = parts[-4]
            method = None
        replicate = parts[-2]

        for record in SeqIO.parse(fname, "fasta"):
            # description actually starts with id
            description = record.description[len(record.id) + 1 :]
            props = dict(pair.split(":") for pair in description.split("|"))

            # extract properties
            freq = props.get("freq")

            if freq is None:
                freq = props.get("Freq")

            # finalize
            tmp.append(
                {
                    "method": method,
                    "params": params,
                    "replicate": replicate,
                    "sequence": str(record.seq),
                    "frequency": float(freq),
                }
            )

    return pd.DataFrame(tmp)


def read_haplostats(haplostats_list):
    df_list = []
    for fname in tqdm(haplostats_list, desc="Read haplostat files"):
        parts = str(fname).split("/")
        params = parts[-4]
        replicate = parts[-2]

        tmp = pd.read_csv(fname)
        tmp["params"] = params
        tmp["replicate"] = replicate

        df_list.append(tmp)

    return pd.concat(df_list)


def read_runstats(runstatus_list):
    tmp = []
    for fname in tqdm(runstatus_list, desc="Read runstatus files"):
        parts = str(fname).split("/")
        params = parts[-5]
        method = parts[-4]
        replicate = parts[-2]

        status = fname.read_text()

        tmp.append(
            {
                "params": params,
                "method": method,
                "replicate": replicate,
                "status": status if len(status) > 0 else "success",
            }
        )

    return pd.DataFrame(tmp)


def overview_plots(df_haplo, dname_out):
    if df_haplo.empty:
        print("Warning: df_haplo is empty")
        return

    df_haplo["seq_len"] = df_haplo["sequence"].str.len()
    df_long = pd.melt(df_haplo, id_vars=["method", "params", "replicate"]).assign(
        params=lambda x: x["params"].str.replace("__", "\n")
    )
    df_long = df_long[df_long["variable"] != "sequence"]

    g = sns.catplot(
        data=df_long,
        x="params",
        y="value",
        hue="method",
        col="variable",
        kind="box",
        sharey=False,
        height=10,
    )
    g.map_dataframe(
        sns.stripplot, x="params", y="value", hue="method", color="k", dodge=True
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
                "--unique-mapping",
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
            parts = str(fname_contigs).split("/")
            if len(parts) == 7:
                _, _, params, method, _, replicate, _ = parts
            elif len(parts) == 8:  # for multi workflow
                _, _, _, params, method, _, replicate, _ = parts
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
            sns.stripplot(data=group, y=col, x="method", color=".25", dodge=True, ax=ax)

            fig.tight_layout()
            fig.savefig(dname_out / f"{col}__{params}.pdf")


def sequence_embedding(df_pred, df_true, dname_out):
    mds_dir = dname_out / "mds_plots"
    mds_dir.mkdir(parents=True)

    for (params, replicate), df_pred_grpd in tqdm(
        df_pred.groupby(["params", "replicate"]), desc="Compute MDS"
    ):
        df_true_grpd = df_true[
            (df_true["params"] == params) & (df_true["replicate"] == replicate)
        ]

        # subsample large results
        max_num = 50

        df_pred_grpd = df_pred_grpd.copy()
        df_pred_grpd = (
            df_pred_grpd.groupby("method")
            .apply(lambda x: x.sample(n=min(len(x), max_num)))
            .reset_index(drop=True)
        )

        # compute dissimilarities
        sequence_list = (
            df_pred_grpd["sequence"].tolist() + df_true_grpd["sequence"].tolist()
        )

        mat = np.zeros(shape=(len(sequence_list), len(sequence_list)))
        for i, seq1 in enumerate(tqdm(sequence_list, leave=False)):
            for j, seq2 in enumerate(tqdm(sequence_list, leave=False)):
                if i >= j:
                    continue

                mat[i, j] = editdistance.eval(seq1, seq2)

        mat = np.triu(mat) + np.tril(mat.T, 1)  # mirror to make symmetric

        # do MDS
        embedding = manifold.MDS(n_components=2, dissimilarity="precomputed")
        mat_trans = embedding.fit_transform(mat)

        df = pd.concat(
            [
                pd.DataFrame(mat_trans, columns=["MDS0", "MDS1"]),
                pd.concat([df_pred_grpd, df_true_grpd], axis=0, ignore_index=True),
            ],
            axis=1,
        )
        df["method"] = df["method"].apply(lambda x: "ground_truth" if x is None else x)

        # plot result
        fig, ax = plt.subplots(figsize=(8, 6))

        sns.scatterplot(data=df, x="MDS0", y="MDS1", hue="method", ax=ax)

        fig.savefig(mds_dir / f"sequence_mds_{params}_{replicate}.pdf")


def compute_pr(df_pred, df_true, thres=0.05):
    @functools.lru_cache(None)
    def compute_dist(seq1, seq2):
        dist = editdistance.eval(seq1, seq2)
        rel = dist / max(len(seq1), len(seq2))
        return rel

    tmp = []
    for (method, params, replicate), df_group in tqdm(
        df_pred.groupby(["method", "params", "replicate"]), desc="Compute PR"
    ):
        tp = 0
        fp = 0
        fn = 0

        # true positive: predicted seq appears in ground truth
        # false positive: predicted seq does not appear in ground truth
        for row in tqdm(df_group.itertuples(), total=df_group.shape[0], leave=False):
            ser_dist = df_true["sequence"].apply(
                lambda x: compute_dist(x, row.sequence)
            )
            passed_thres = (ser_dist <= thres).any()

            if passed_thres:
                tp += 1
            else:
                fp += 1

        # false negative: ground truth sequence was not predicted
        # single prediction should not map to multiple ground truth seqs
        df_cur = df_group.copy()
        for row in tqdm(df_true.itertuples(), total=df_true.shape[0], leave=False):
            ser_dist = df_cur["sequence"].apply(lambda x: compute_dist(x, row.sequence))
            passed_thres = (ser_dist <= thres).any()

            if not passed_thres:
                fn += 1
            else:
                # remove current prediction
                df_cur = df_cur.drop(ser_dist.idxmin())

        # finalize
        tmp.append(
            {
                "method": method,
                "params": params,
                "replicate": replicate,
                "tp": tp,
                "fp": fp,
                "fn": fn,
                "precision": tp / (tp + fp),
                "recall": tp / (tp + fn),
            }
        )

    # set column dtypes
    df_pr = pd.DataFrame(tmp)
    df_pr["method"] = pd.Categorical(
        df_pr["method"], categories=sorted(snakemake.params.method_list_global)
    )

    return df_pr


def plot_pr(df_pr, df_stats, dname_out):
    # prepare data
    diversity_column_list = ["population_nucleotide_diversity", "mean_position_shannon"]

    df_m = df_pr.merge(df_stats, on=["params", "replicate"]).assign(
        params=lambda x: x["params"].str.replace("__", "\n")
    )

    # helper functions
    def do_plot(df, x, y, fname):
        fig, ax = plt.subplots()

        sns.boxplot(data=df, x=x, y=y, hue="method", ax=ax)
        sns.swarmplot(
            data=df,
            x=x,
            y=y,
            hue="method",
            dodge=True,
            clip_on=False,
            linewidth=1,
            edgecolor="gray",
            ax=ax,
        )

        ax.set_ylim(0, 1)

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[:2], labels[:2])

        fig.tight_layout()
        fig.savefig(fname)

    # plots
    do_plot(df_m, "params", "precision", dname_out / "overview_precision.pdf")
    do_plot(df_m, "params", "recall", dname_out / "overview_recall.pdf")

    for diversity_column in diversity_column_list:
        do_plot(
            df_m,
            diversity_column,
            "precision",
            dname_out / f"overview_precision_{diversity_column}.pdf",
        )
        do_plot(
            df_m,
            diversity_column,
            "recall",
            dname_out / f"overview_recall_{diversity_column}.pdf",
        )


def main(
    predicted_haplos_list, true_haplos_list, haplostats_list, runstatus_list, dname_out
):
    dname_out.mkdir(parents=True)

    # read data
    df_pred = read_fasta_files(predicted_haplos_list)
    df_true = read_fasta_files(true_haplos_list, with_method=False)
    df_true["method"] = "ground_truth"

    df_stats = read_haplostats(haplostats_list)
    df_runstats = read_runstats(runstatus_list)

    # quick stats
    print("Run status")
    print(df_runstats.groupby("method")["status"].value_counts())

    print("Haplotype counts per method")
    print(df_pred["method"].value_counts())

    # create plots
    overview_plots(df_pred, dname_out)

    # precision/recall
    df_pr = compute_pr(df_pred, df_true)
    plot_pr(df_pr, df_stats, dname_out)

    # quast stuff
    df_quast = run_metaquast(
        predicted_haplos_list, true_haplos_list, dname_out / "quast"
    )
    plot_quast(df_quast, dname_out)

    # MDS
    sequence_embedding(df_pred, df_true, dname_out)

    # save raw results
    csv_dir = dname_out / "csv_files"
    csv_dir.mkdir(parents=True, exist_ok=True)

    df_pred.to_csv(csv_dir / "predictions.csv.gz")
    df_true.to_csv(csv_dir / "ground_truth.csv.gz")
    df_stats.to_csv(csv_dir / "data_stats.csv")
    df_runstats.to_csv(csv_dir / "run_stats.csv")
    df_pr.to_csv(csv_dir / "pr_results.csv")
    df_quast.to_csv(csv_dir / "quast_results.csv")


if __name__ == "__main__":
    main(
        [Path(e) for e in snakemake.input.predicted_haplos_list],
        [Path(e) for e in snakemake.input.true_haplos_list],
        [Path(e) for e in snakemake.input.haplostats_list],
        [Path(e) for e in snakemake.input.runstatus_list],
        Path(snakemake.output.dname_out),
    )
