#!/usr/bin/env python3
"""
Computation of various diversity indices for the underlying sample following the
review: https://doi.org/10.1016/j.coviro.2021.06.002

"""
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
    n_reads_with_mut = df_temp["frequency"] * df_temp["coverage"]
    n_reads_with_mut = df_temp["rvar"] + df_temp["fvar"]
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
        if N == 0:
            continue
        freq = df_temp["frequency"].to_numpy()
        ref_freq = 1 - freq.sum()

        position_pnd = freq**2
        postion_pi = (1 - (position_pnd.sum() + ref_freq**2)) * N / (N - 1)

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


def convert_vcf(fname):
    """Convert VCF to JSON."""
    output = []

    if os.path.getsize(fname) == 0:
        print(f'Empty VCF: "{fname}"')
        return output

    print(f'Parsing VCF: "{fname}"')
    with open(fname) as fd:
        vcf_reader = vcf.Reader(fd)

        # check caller
        caller_source = vcf_reader.metadata["source"][0].lower()
        if caller_source.startswith("lofreq"):
            mode = "lofreq"
        elif caller_source.startswith("shorah"):
            mode = "shorah"
        else:
            raise RuntimeError(f"Invalid variant caller: {caller_source}")

        # output dataframe
        cols_snv = [
            "position",
            "reference",
            "variant",
            "frequency",
            "tvar",
            "fvar",
            "rvar",
            "ftot",
            "rtot",
            "coverage",
        ]
        df_snv = pd.DataFrame(columns=cols_snv)

        # parse records
        for record in vcf_reader:
            if mode == "lofreq":
                freq = round(record.INFO["AF"], 3)
                fvar = record.INFO["DP4"][2]  # counts variant forward reads
                rvar = record.INFO["DP4"][3]  # counts variant reverse reads
                rtot = (
                    record.INFO["DP4"][1] + record.INFO["DP4"][3]
                )  # counts reverse reads
                ftot = (
                    record.INFO["DP4"][0] + record.INFO["DP4"][2]
                )  # counts forward reads

            elif mode == "shorah":
                freq = round(
                    np.mean(
                        [v for k, v in record.INFO.items() if k.startswith("Freq")]
                    ),
                    3,
                )
                fvar = record.INFO["Fvar"]  # counts variant forward reads
                rvar = record.INFO["Rvar"]  # counts variant reverse reads
                rtot = record.INFO["Rtot"]  # counts reference reverse reads
                ftot = record.INFO["Ftot"]  # counts reference forward reads

            df_snv = df_snv.append(
                {
                    "position": record.POS,
                    "reference": record.REF,
                    "variant": [v.sequence for v in record.ALT],
                    "frequency": freq,
                    "tvar": fvar + rvar,
                    "fvar": fvar,
                    "rvar": rvar,
                    "ftot": ftot,
                    "rtot": rtot,
                    "coverage": ftot + rtot,
                },
                ignore_index=True,
            )

        id = record.CHROM
    return df_snv, id


def load_reference_seq(reference_file):
    for seq in skbio.io.read(reference_file, format="fasta"):
        return seq


def main(fname_snv_in, fname_reference, output_diversity_csv, output_shannon_csv):
    """
    Compute various diversity indices for underlying sample.
    Writes a csv-file with computations.

    Parameters
    ----------
    fname_snv_in:
        Absolute path to snv.vcf created by lofreq or ShoRAH.
    ref_seq_length:
        Length of the reference sequence.
    output_diversity_csv:
        Absolute path to diversity-csv.
    output_shannon_csv:
        Absolute path to position-wise Shannon-entropy csv.

    Output
    ------
    diversity_measures.csv:
        Listing all computed diversity measures.
    position_shannon_entropy.csv:
        Position-wise Shannon entropy, if non-zero.
    """
    # get length of reference sequence
    ref_seq_length = len(load_reference_seq(fname_reference))

    # Parse snv.vcf
    df_snv, id = convert_vcf(fname_snv_in)

    # Sample, patient, date information
    general_info = {
        "sample": fname_snv_in.split("/variants")[0].split("/")[-3],
        "patient": fname_snv_in.split("/variants")[0].split("/")[-2],
        "date": fname_snv_in.split("/variants")[0].split("/")[-1],
    }

    # prepare dict collecting all the diversity measures
    out_dict = {"id": id, "length": ref_seq_length}

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
            "population_nucleotide_diverstiy": population_nucleotide_diversity(
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

    out_dict.update(general_info)

    # save to dataframe
    df_diveristy = pd.DataFrame(columns=list(out_dict.keys()))
    df_diveristy = df_diveristy.append(out_dict, ignore_index=True)
    df_diveristy.to_csv(output_diversity_csv, index=False)

    # position-wise Shannon entropy
    pos_shannon_dict = general_info
    for i in range(ref_seq_length):
        if position_Shannon_entropy(df_snv, i) != 0:
            pos_shannon_dict.update({i: position_Shannon_entropy(df_snv, i)})
    df_pos_shannon = pd.DataFrame(columns=list(pos_shannon_dict.keys()))
    df_pos_shannon = df_pos_shannon.append(pos_shannon_dict, ignore_index=True)
    df_pos_shannon.to_csv(output_shannon_csv, index=False)


if __name__ == "__main__":
    main(
        snakemake.input.fnames_samples_snvs_vcf,
        snakemake.input.fname_ref,
        snakemake.output.diversity_csv,
        snakemake.output.shannon_csv,
    )
