from pathlib import Path

import humanize
import pandas as pd
from cyvcf2 import VCF
import os
import seaborn as sns
from matplotlib.ticker import FuncFormatter


@FuncFormatter
def duration_fmt(x, pos):
    return humanize.naturaldelta(x)


def convert_vcf(df):
    # TODO: check what happens with indels
    variant_list = set()

    for idx_row, variant in df.iterrows():
        for base in variant['Alt']:
            zero_based_pos = int(variant['Pos']) - 1  # VCF is 1-based
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


def performance_plots(vcf_list, groundtruth_list, posterior_threshold):
    # compute performance
    colnames=['Chromosome', 'Pos', 'Ref', 'Alt', 'Frq', 'Pst', 'Rvar', 'Fvar', 'Rtot', 'Ftot', 'Qval']
    tmp = []
    fps_tmp = []
    for fname_vcf, fname_groundtruth in zip(vcf_list, groundtruth_list):

        # we want the file SNVs_0.010000.tsv as the posterior filter was not applied here
        fname_SNVs_correct = str(fname_vcf).split("snvs_.vcf")[0] + "work/snv/SNVs_0.010000.tsv"
        if os.path.isfile(fname_SNVs_correct):
            df_predicted = pd.read_csv(fname_SNVs_correct, names=colnames, header=None ,sep="\t")
        else:
            fname_SNVs_correct = str(fname_vcf).split("snvs_.vcf")[0] + "work/snv/SNVs_0.010000_final.csv"
            df_predicted = pd.read_csv(fname_SNVs_correct)
            df_predicted = df_predicted.rename(columns={'Var': 'Alt',
                                                        'Pst2': 'Pst'})

        parts = str(fname_vcf).split("/")

        if len(parts) == 7:
            _, _, params, method, _, replicate, _ = parts
        elif len(parts) == 8:  # for multi workflow
            _, _, _, params, method, _, replicate, _ = parts


        #filter dataframe according to posterior_threshold
        df_predicted['Pst'] = pd.to_numeric(df_predicted['Pst'], errors='coerce')
        df_predicted = df_predicted.dropna(subset=['Pst'])
        df_predicted = df_predicted[df_predicted['Pst']>posterior_threshold]

        true_variants = convert_groundtruth(fname_groundtruth)
        predicted_variants = convert_vcf(df_predicted)

        if len(true_variants) == 0:
            # no true variants
            # Goal: Count the false positives
            fp = len(predicted_variants)
            fps_tmp.append(
                {
                    "fname_vcf": fname_vcf,
                    "method": method,
                    "params": params,
                    "replicate": replicate,
                    "fp": fp,
                }
            )

        precision, recall, f1 = compute_performance(true_variants, predicted_variants)

        tmp.append(
            {
                "method": method,
                "params": params,
                "replicate": replicate,
                "posterior_threshold": posterior_threshold,
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
    return df_perf



def main(vcf_list, groundtruth_list, fname_out):

    posterior_thresholds = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
    dfs_tmp = []

    for posterior_threshold in posterior_thresholds:

        df_pos_thres = performance_plots(
                                vcf_list,
                                groundtruth_list,
                                posterior_threshold
                                )

        dfs_tmp.append(df_pos_thres)

    pd.concat(dfs_tmp).to_csv(fname_out)

if __name__ == "__main__":
    main(
        [Path(e) for e in snakemake.input.vcf_list],
        [Path(e) for e in snakemake.input.groundtruth_list],
        snakemake.output.fname_out,
    )
