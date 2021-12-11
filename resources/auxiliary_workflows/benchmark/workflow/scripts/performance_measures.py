from pathlib import Path

import pandas as pd
from cyvcf2 import VCF

import seaborn as sns


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
    f1 = 2 * (precision * recall) / (precision + recall)

    return precision, recall, f1


def main(vcf_list, groundtruth_list, dname_out):
    dname_out.mkdir(parents=True)

    # compute performance
    tmp = []
    for fname_vcf, fname_groundtruth in zip(vcf_list, groundtruth_list):
        _, _, params, method, _, replicate, _ = str(fname_vcf).split("/")

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

    g = sns.catplot(
        data=df_long,
        x="params",
        y="value",
        hue="method",
        row="variable",
        kind="box",
    )
    g.set(ylim=(0, 1))
    g.savefig(dname_out / "overview.pdf")


if __name__ == "__main__":
    main(
        [Path(e) for e in snakemake.input.vcf_list],
        [Path(e) for e in snakemake.input.groundtruth_list],
        Path(snakemake.output.dname_out),
    )
