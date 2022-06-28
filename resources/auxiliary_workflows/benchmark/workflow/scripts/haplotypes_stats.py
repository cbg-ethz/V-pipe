#!/usr/bin/env python3

import sys
import os
import vcf
import pandas as pd
import numpy as np
import skbio
from scipy.stats import sem


def number_of_polymorphisms(df_mutations, minor_allele_frequency=0):
    df_temp = df_mutations[df_mutations["frequency"] >= minor_allele_frequency]
    variant_positions = df_temp["position"].unique()
    return len(variant_positions)


def list_polymorphic_sites(df_mutations, minor_allele_frequency=0):
    df_temp = df_mutations[df_mutations["frequency"] >= minor_allele_frequency]
    variant_positions = df_temp["position"].unique()
    return variant_positions


def number_of_mutations(df_mutations, minimum_frequency=0):
    df_temp = df_mutations[df_mutations["frequency"] >= minimum_frequency]
    n_reads_with_mut = df_temp["frequency"] * pd.to_numeric(
        df_temp["coverage"], downcast="float"
    )
    n_reads_with_mut = df_temp["tvar"]
    return n_reads_with_mut.sum()


def mutation_spectrum(df_mutations, bins):
    counts = np.zeros(len(bins))
    for i in range(len(bins) - 1):
        df_temp = df_mutations[df_mutations["frequency"] > bins[i]]
        df_temp = df_temp[df_temp["frequency"] <= bins[i + 1]]
        n_reads_with_mut = df_temp["frequency"] * df_temp["coverage"]
        counts[i + 1] = n_reads_with_mut.sum()
    return list(counts)


def population_nucleotide_diversity(df_mutations, length):
    # only the positions with mutations are needed
    pi = 0
    for position_temp in list_polymorphic_sites(df_mutations, minor_allele_frequency=0):
        df_temp = df_mutations[df_mutations["position"] == position_temp]
        N = df_temp["coverage"].unique()[0]
        position_pnd = 0
        pi = df_temp["tvar"] * (df_temp["tvar"] - 1)
        postion_pi = (N * (N - 1) - pi.sum()) / (N * (N - 1))

        pi += postion_pi

    return float(pi / length)


def position_Shannon_entropy(df_mutations, position):
    df_temp = df_mutations[df_mutations["position"] == position]
    position_shannon = 0
    sum_fraction = 0
    var_shannon = df_temp["frequency"].apply(lambda x: x * np.log(x) if x > 0 else 0)

    sum_fraction = df_temp["frequency"].sum()
    position_shannon = var_shannon.sum()

    # add the reference base summand
    if 1 - sum_fraction > 0:
        position_shannon += (1 - sum_fraction) * np.log(1 - sum_fraction)

    return -position_shannon


def mean_pos_Shannon_entropy(df_mutations, length):
    entropy = 0
    for position_temp in list_polymorphic_sites(df_mutations, minor_allele_frequency=0):
        entropy += position_Shannon_entropy(df_mutations, position_temp)
    return entropy / length


def load_reference_seq(reference_file):
    for seq in skbio.io.read(reference_file, format="fasta"):
        return seq


def main(fname_snv_in, fname_reference, coverage, fname_out):
    """
    Compute various diversity indices.
    """
    # get length of reference sequence
    ref_seq_length = len(load_reference_seq(fname_reference))

    # Parse ground_truth
    df_snv = pd.read_csv(fname_snv_in)
    df_snv["coverage"] = float(coverage)
    df_snv["frequency"] = pd.to_numeric(df_snv["frequency"], downcast="float")
    df_snv["tvar"] = df_snv["frequency"].apply(lambda x: int(x * float(coverage)))

    # prepare dict collecting all the diversity measures
    out_dict = {"genome_length": ref_seq_length}

    # number of mutations with different minor allele frequency
    out_dict.update(
        {
            "n_mutations_minFrq_0": number_of_mutations(df_snv, minimum_frequency=0),
            "n_mutations_minFrq_1": number_of_mutations(df_snv, minimum_frequency=0.01),
            "n_mutations_minFrq_5": number_of_mutations(df_snv, minimum_frequency=0.05),
        }
    )

    # sum mutation frequencies
    out_dict.update({"sum_mutation_frq": df_snv["frequency"].sum()})

    # mean mutation frequencies
    out_dict.update({"mean_mutation_frq": df_snv["frequency"].mean()})

    # the standard error of the mean (SEM) mutation frequency
    out_dict.update({"sem_mutation_frq": sem(df_snv["frequency"].to_numpy())})

    # population nucleotide diversity
    out_dict.update(
        {
            "population_nucleotide_diversity": population_nucleotide_diversity(
                df_snv, ref_seq_length
            )
        }
    )

    # mean position-wise Shannon entropy
    out_dict.update(
        {"mean_position_shannon": mean_pos_Shannon_entropy(df_snv, ref_seq_length)}
    )

    # mutation spectrum
    bins = np.arange(0, 1, 0.05)
    out_dict.update({"mutation_spectrum_bins": list(bins)})
    out_dict.update({"mutation_spectrum": mutation_spectrum(df_snv, bins)})

    # save to dataframe
    df_diveristy = pd.DataFrame(columns=list(out_dict.keys()))
    df_diveristy = df_diveristy.append(out_dict, ignore_index=True)
    df_diveristy.to_csv(fname_out, index=False)


if __name__ == "__main__":
    main(
        snakemake.input.fname_groundtruth,
        snakemake.input.fname_reference,
        snakemake.wildcards.coverage,
        snakemake.output.fname,
    )
