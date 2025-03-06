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


def convert_vcf(fname):
    # TODO: check what happens with indels
    variant_list = set()
    try:
        vcf = VCF(fname)
    except OSError:
        return set()

    for variant in VCF(fname):
        for base in variant.ALT:
            zero_based_pos = variant.POS - 1  # VCF is 1-based
            variant_list.add(f"{zero_based_pos}{base}")
    return variant_list


def convert_groundtruth(fname):
    df = pd.read_csv(fname, index_col=0)
    return set((df["position"].astype(str) + df["variant"]).tolist())


def mutation_calls_details(vcf_list, groundtruth_list):
    # compute performance
    tmp = []
    for fname_vcf, fname_groundtruth in zip(vcf_list, groundtruth_list):
        parts = str(fname_vcf).split("/")

        if len(parts) == 7:
            _, _, params, method, _, replicate, _ = parts
        elif len(parts) == 8:  # for multi workflow
            _, _, _, params, method, _, replicate, _ = parts

        predicted_variants = convert_vcf(fname_vcf)

        # iter through ground truth mutations
        for iter_row, gt_row in pd.read_csv(fname_groundtruth, index_col=0).iterrows():
            true_variant = str(gt_row["position"]) + gt_row["variant"]
            mutation_type = str(gt_row["type"])
            if method.startswith(("lofreq", "shorah")):
                is_false_negative = True
                for variant in VCF(fname_vcf):
                    for idx_base, base in enumerate(variant.ALT):
                        zero_based_pos = variant.POS - 1  # VCF is 1-based
                        predicted_variant = f"{zero_based_pos}{base}"

                        if true_variant == predicted_variant:
                            is_false_negative = False
                            # get pred_freq
                            if method.startswith("lofreq"):
                                posterior = 1.0
                                freq = variant.INFO.get("AF")  # vi.split(',')[idx_base]
                            elif method.startswith("shorah"):
                                posterior = 1.0
                                freq = (
                                    int(variant.INFO.get("Fvar"))
                                    + int(variant.INFO.get("Rvar"))
                                ) / (
                                    int(variant.INFO.get("Ftot"))
                                    + int(variant.INFO.get("Rtot"))
                                )

                            tmp.append(
                                {
                                    "method": method,
                                    "params": params,
                                    "replicate": replicate,
                                    "mutation_type": mutation_type,
                                    "true_positive": 1,
                                    "false_negative": 0,
                                    "haplotype": gt_row["haplotype"],
                                    "frequency_groundtruth": gt_row["frequency"],
                                    "frequency_predicted": freq,
                                    "posterior": posterior,
                                    "position": gt_row["position"],
                                    "variant": gt_row["variant"],
                                }
                            )
                if is_false_negative:
                    tmp.append(
                        {
                            "method": method,
                            "params": params,
                            "mutation_type": mutation_type,
                            "replicate": replicate,
                            "true_positive": 0,
                            "false_negative": 1,
                            "haplotype": gt_row["haplotype"],
                            "frequency_groundtruth": gt_row["frequency"],
                            "frequency_predicted": 0.0,
                            "posterior": 1.0,
                            "position": gt_row["position"],
                            "variant": gt_row["variant"],
                        }
                    )

            elif method.startswith("viloca"):

                colnames = [
                    "Chromosome",
                    "Pos",
                    "Ref",
                    "Alt",
                    "Frq",
                    "Pst",
                    "Rvar",
                    "Fvar",
                    "Rtot",
                    "Ftot",
                    "Qval",
                ]
                # we want the file SNVs_0.010000.tsv as the posterior filter was not applied here
                fname_SNVs_correct = (
                    str(fname_vcf).split("snvs_.vcf")[0] + "work/snv/SNVs_0.010000.tsv"
                )
                df_predicted = pd.read_csv(
                    fname_SNVs_correct, names=colnames, header=None, sep="\t"
                )
                df_predicted["Pst"] = pd.to_numeric(
                    df_predicted["Pst"], errors="coerce"
                )

                for iter_row_pred, pred_row in df_predicted.iterrows():
                    zero_based_pos = int(pred_row["Pos"]) - 1  # VCF is 1-based
                    is_false_negative = True
                    predicted_variant = f"{zero_based_pos}{pred_row['Alt']}"

                    if true_variant == predicted_variant:
                        is_false_negative = False
                        posterior = pred_row["Pst"]
                        freq = pred_row["Frq"]

                        tmp.append(
                            {
                                "method": method,
                                "params": params,
                                "replicate": replicate,
                                "true_positive": 1,
                                "false_negative": 0,
                                "haplotype": gt_row["haplotype"],
                                "frequency_groundtruth": gt_row["frequency"],
                                "frequency_predicted": freq,
                                "posterior": posterior,
                                "position": gt_row["position"],
                                "variant": gt_row["variant"],
                            }
                        )
                if is_false_negative:
                    tmp.append(
                        {
                            "method": method,
                            "params": params,
                            "mutation_type": mutation_type,
                            "replicate": replicate,
                            "true_positive": 0,
                            "false_negative": 1,
                            "haplotype": gt_row["haplotype"],
                            "frequency_groundtruth": gt_row["frequency"],
                            "frequency_predicted": 0.0,
                            "posterior": 1.0,
                            "position": gt_row["position"],
                            "variant": gt_row["variant"],
                        }
                    )

    return pd.DataFrame(tmp)


def main(vcf_list, groundtruth_list, fname_out):

    df = mutation_calls_details(vcf_list, groundtruth_list)
    df.to_csv(fname_out)


if __name__ == "__main__":
    main(
        [Path(e) for e in snakemake.input.vcf_list],
        [Path(e) for e in snakemake.input.groundtruth_list],
        snakemake.output.fname_out,
    )
