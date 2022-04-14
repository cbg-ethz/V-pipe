from pathlib import Path

import pandas as pd
from Bio import SeqIO

import seaborn as sns


def read_fasta_files(fasta_files):
    tmp = []
    for fname in fasta_files:
        _, _, params, method, _, replicate, _ = str(fname).split("/")

        for record in SeqIO.parse(fname, "fasta"):
            tmp.append(
                {
                    "method": method,
                    "params": params,
                    "replicate": replicate,
                    "sequence": str(record.seq),
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


def main(fasta_list, dname_out):
    dname_out.mkdir(parents=True)

    df_haplo = read_fasta_files(fasta_list)
    overview_plots(df_haplo, dname_out)


if __name__ == "__main__":
    main(
        [Path(e) for e in snakemake.input.fasta_list],
        Path(snakemake.output.dname_out),
    )