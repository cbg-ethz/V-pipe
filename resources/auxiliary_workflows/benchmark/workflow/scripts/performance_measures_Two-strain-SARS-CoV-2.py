import pandas as pd
from cyvcf2 import VCF

def convert_vcf(fname):
    # TODO: check what happens with indels
    variant_list = set()
    tmp = []
    for variant in VCF(fname):
        for base in variant.ALT:
            one_based_pos = variant.POS  # VCF is 1-based
            freq = variant.AF
            variant_list.add(f"{one_based_pos}{base}")
            tmp.append(
                "one_based_pos": one_based_pos,
                "freq": freq,
                "base": base,
            )
    df = pd.DataFrame(tmp)
    return variant_list, df

def convert_groundtruth(fname):
    df = pd.read_csv(fname, index_col=0)
    return set((df["position"].astype(str) + df["variant"]).tolist())

def performance_evaluation(vcf_list, groundtruth_list, dname_out):

    #TODO: how to make sure that we are also catching the deletions
    # probably needs to be reported differnetlich in the ground_truth such that there is a match
    tmp = []
    true_variants = convert_groundtruth(groundtruth_list[0])

    for variant in true_variants:
        for fname_vcf in vcf_list:
            parts = str(fname_vcf).split("/")

            if len(parts) == 7:
                _, _, params, method, _, replicate, _ = parts
            elif len(parts) == 8: # for multi workflow
                _, _, _, params, method, _, replicate, _ = parts

            predicted_variants, df_predicted = convert_vcf(fname_vcf)

            position = int(variant[:-1])
            base = str(variant[-1])

            if variant in predicted_variants:

                df = df_predicted[df_predicted['one_based_pos'] == position]
                df = df[df['base']== base]
                freq = df['freq'].values[0]

            else:
                freq = 0


            tmp.append(
                {
                    "position": position,
                    "variant": base,
                    "params": params,
                    "replicate": replicate,
                    "freq": freq,
                }
            )

    df_perf = pd.DataFrame(tmp)
    df_long = pd.melt(df_perf, id_vars=["method", "params", "replicate"]).assign(
        params=lambda x: x["params"].str.replace("_", "\n")
    )
    df_long.to_csv(fname_out)


def main(vcf_list, groundtruth_list, fname_out):

    dname_out.mkdir(parents=True)
    performance_evaluation(vcf_list, groundtruth_list, fname_out)


if __name__ == "__main__":
    main(
        [Path(e) for e in snakemake.input.vcf_list],
        [Path(e) for e in snakemake.input.groundtruth_list],
        Path(snakemake.output.fname_out),
    )
