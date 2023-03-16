from pathlib import Path

import humanize
import pandas as pd
from cyvcf2 import VCF

import seaborn as sns
from matplotlib.ticker import FuncFormatter


@FuncFormatter
def duration_fmt(x, pos):
    return humanize.naturaldelta(x)


def convert_vcf(fname):
    # TODO: check what happens with indels
    variant_list = set()
    for variant in VCF(fname):
        for base in variant.ALT:
            zero_based_pos = variant.POS - 1  # VCF is 1-based
            variant_list.add(f"{zero_based_pos}{base}")
    return variant_list


def convert_groundtruth(fname):
    df = pd.read_csv(fname, index_col=0)
    return set((df["position"].astype(str) + df["variant"]).tolist())


def compute_performance(true_variants, predicted_variants):
    # no true variants and also no predicted variants
    if len(true_variants) == 0 and len(predicted_variants) == 0:
        tp, fp, fn = 0, 0, 0
        precision = 1
        recall = 1
        f1 = 2 * (precision * recall) / (precision + recall)
    elif len(true_variants) == 0 and len(predicted_variants) > 0:
        tp, fn = 0, 0
        fp = len(predicted_variants)
        precision = tp / (tp + fp)
        recall = 0  # 0 / 0
        f1 = 0
    elif len(true_variants) > 0 and len(predicted_variants) == 0:
        tp, fp = 0, 0
        fn = len(true_variants)
        precision = 0  # 0 / 0
        recall = tp / (tp + fn)
        f1 = 0
    else:
        # count true/false positives/negatives
        tp, fp, fn = 0, 0, 0
        for variant in true_variants | predicted_variants:
            if variant in true_variants and variant in predicted_variants:
                tp += 1
            elif variant in true_variants and variant not in predicted_variants:
                fn += 1
            elif variant not in true_variants and variant in predicted_variants:
                fp += 1
            elif variant not in true_variants and variant not in predicted_variants:
                # tn += 1
                pass
            else:
                raise RuntimeError("Woopsie")

        # compute performances
        precision = tp / (tp + fp)
        recall = tp / (tp + fn)
        if tp == 0:
            f1 = 0
        else:
            f1 = 2 * (precision * recall) / (precision + recall)

    return precision, recall, f1


def performance_plots(vcf_list, groundtruth_list, dname_out):
    # compute performance
    tmp = []
    for fname_vcf, fname_groundtruth in zip(vcf_list, groundtruth_list):
        parts = str(fname_vcf).split("/")

        if len(parts) == 7:
            _, _, params, method, _, replicate, _ = parts
        elif len(parts) == 8:  # for multi workflow
            _, _, _, params, method, _, replicate, _ = parts

        true_variants = convert_groundtruth(fname_groundtruth)
        predicted_variants = convert_vcf(fname_vcf)

        precision, recall, f1 = compute_performance(true_variants, predicted_variants)

        tmp.append(
            {
                "method": method,
                "params": params,
                "replicate": replicate,
                "precision": precision,
                "recall": recall,
                "f1": f1,
            }
        )
    df_perf = pd.DataFrame(tmp)

    # plot overview
    df_long = pd.melt(df_perf, id_vars=["method", "params", "replicate"]).assign(
        params=lambda x: x["params"].str.replace("_", "\n")
    )
    df_long.to_csv(dname_out / "performance.csv")

    g = sns.catplot(
        data=df_long,
        x="params",
        y="value",
        hue="method",
        col="variable",
        kind="box",
    )
    g.set(ylim=(0, 1))
    g.savefig(dname_out / "performance_boxplot.pdf")


def runtime_plots(benchmark_list, dname_out):
    # gather benchmark information
    tmp = []
    for fname in benchmark_list:
        parts = str(fname).split("/")
        if len(parts) == 7:
            _, _, params, method, _, replicate, _ = parts
        elif len(parts) == 8:  # for multi workflow
            _, _, _, params, method, _, replicate, _ = parts

        df_tmp = pd.read_csv(fname, sep="\t")
        df_tmp["method"] = method
        df_tmp["params"] = params
        df_tmp["replicate"] = replicate

        tmp.append(df_tmp)
    df_bench = pd.concat(tmp).replace("-", pd.NA)

    # plot
    df_long = (
        pd.melt(df_bench, id_vars=["method", "params", "replicate"])
        .assign(params=lambda x: x["params"].str.replace("_", "\n"))
        .query("variable == 's'")
    )
    df_long.to_csv(dname_out / "runtime.csv")

    g = sns.catplot(
        data=df_long,
        x="params",
        y="value",
        hue="method",
        col="variable",
        kind="box",
    )

    for ax in g.axes.flat:
        ax.yaxis.set_major_formatter(duration_fmt)

    g.savefig(dname_out / "runtime_boxplot.pdf")


def main(vcf_list, groundtruth_list, benchmark_list, dname_out):
    dname_out.mkdir(parents=True)

    performance_plots(vcf_list, groundtruth_list, dname_out)
    runtime_plots(benchmark_list, dname_out)


if __name__ == "__main__":
    main(
        [Path(e) for e in snakemake.input.vcf_list],
        [Path(e) for e in snakemake.input.groundtruth_list],
        [Path(e) for e in snakemake.input.benchmark_list],
        Path(snakemake.output.dname_out),
    )
