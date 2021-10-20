#!/usr/bin/env python3
"""
Computation of various diversity indices for the underlying sample following the
review: https://doi.org/10.1016/j.coviro.2021.06.002

TODO: Option to only include positions with coverage>1000
"""
import sys
import os
import vcf
import pandas as pd
import numpy as np
import skbio
from scipy.stats import sem


def number_of_polymorphisms(df_mutations, minor_allele_frequency=0):
    df_temp = df_mutations[df_mutations['frequency']>=minor_allele_frequency]
    variant_positions = df_temp['position'].unique()
    return len(variant_positions)

def list_polymorphic_sites(df_mutations, minor_allele_frequency=0):
    df_temp = df_mutations[df_mutations['frequency']>=minor_allele_frequency]
    variant_positions = df_temp['position'].unique()
    return variant_positions

def number_of_mutations(df_mutations, minimum_frequency=0):
    df_temp = df_mutations[df_mutations['frequency']>=minimum_frequency]
    return df_temp.shape[0]

def population_nucleotide_diversity(df_mutations, length):
    # only the positions with mutations are needed
    pi = 0
    for position_temp in list_polymorphic_sites(df_mutations, minor_allele_frequency=0):
        df_temp = df_mutations[df_mutations['position']== position_temp]
        N = df_temp['coverage'].unique()[0]
        position_pnd =0
        df_temp['pi']= df_temp['tvar']*(df_temp['tvar']-1)
        postion_pi = (N*(N-1)-df_temp['pi'].sum())/(N*(N-1))

        pi+=postion_pi

    return pi/length

def position_Shannon_entropy(df_mutations,position):
    df_temp= df_mutations[df_mutations['position']==position]
    position_shannon = 0
    sum_fraction = 0
    df_temp['var_shannon']= df_temp['frequency'].apply(lambda x: x*np.log(x))

    sum_fraction = df_temp['frequency'].sum()
    position_shannon = df_temp['var_shannon'].sum()

    # add the reference base summand
    position_shannon+=(1-sum_fraction)*np.log(1-sum_fraction)

    return - position_shannon

def mean_pos_Shannon_entropy(df_mutations, length):
    entropy = 0
    for position_temp in list_polymorphic_sites(df_mutations, minor_allele_frequency=0):
        entropy+=position_Shannon_entropy(df_mutations,position_temp)
    return entropy/length

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
        caller_source = vcf_reader.metadata['source'][0].lower()
        if caller_source.startswith('lofreq'):
            mode = 'lofreq'
        elif caller_source.startswith('shorah'):
            mode = 'shorah'
        else:
            raise RuntimeError(f'Invalid variant caller: {caller_source}')

        # output dataframe
        cols_snv = ["position","reference","variant","frequency","tvar","fvar",
                    "rvar","ftot","rtot","coverage"]
        df_snv = pd.DataFrame(columns=cols_snv)

        # parse records
        for record in vcf_reader:
            if mode == 'lofreq':
                freq = round(record.INFO['AF'], 3)
                fvar = record.INFO['DP4'][2] # counts variant forward reads
                rvar = record.INFO['DP4'][3] # counts variant reverse reads
                rtot = record.INFO['DP4'][1] + record.INFO['DP4'][3] # counts reverse reads
                ftot = record.INFO['DP4'][0] + record.INFO['DP4'][2] # counts forward reads

            elif mode == 'shorah':
                freq = round(np.mean(
                    [v for k, v in record.INFO.items() if k.startswith("Freq")]
                ), 3)
                fvar = record.INFO['Fvar'] # counts variant forward reads
                rvar = record.INFO['Rvar'] # counts variant reverse reads
                rtot = record.INFO['Rtot'] # counts reference reverse reads
                ftot = record.INFO['Ftot'] # counts reference forward reads

            df_snv=df_snv.append(
                {
                    "position": record.POS,
                    "reference": record.REF,
                    "variant": [v.sequence for v in record.ALT],
                    "frequency": freq,
                    "tvar": fvar+rvar,
                    "fvar": fvar,
                    "rvar": rvar,
                    "ftot": ftot,
                    "rtot": rtot,
                    "coverage": ftot+rtot
                },
                ignore_index=True
            )

        id = record.CHROM
    return df_snv, id

def load_reference_seq(reference_file):
    for seq in skbio.io.read(reference_file, format='fasta'):
         return seq

def main(fname_snv_in, fname_reference ,output_dir):
    """
    Compute various diversity indices for underlying sample.
    Writes a csv-file with computations.

    Parameters
    ----------
    fname_snv_in:
        Absolute path to snv.vcf created by lofreq or ShoRAH.
    ref_seq_length:
        Length of the reference sequence.
    output_dir:
        Path to output directory.

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

    # prepare dict collecting all the diversity measures
    out_dict = {"id": id, "length": ref_seq_length}

    # number of mutations with different minor allele frequency
    out_dict.update({'n_mutations_minFrq_0': number_of_mutations(df_snv, minimum_frequency =0),
                    'n_mutations_minFrq_1': number_of_mutations(df_snv, minimum_frequency =0.01),
                    'n_mutations_minFrq_5': number_of_mutations(df_snv, minimum_frequency =0.05)})

    # sum mutation frequencies
    out_dict.update({'sum_mutation_frq': df_snv['frequency'].sum()})

    # mean mutation frequencies
    out_dict.update({'mean_mutation_frq': df_snv['frequency'].mean()})

    # the standard error of the mean (SEM) mutation frequency
    out_dict.update({'sem_mutation_frq': sem(df_snv['frequency'].to_numpy())})

    # population nucleotide diversity
    out_dict.update({'n_population_nucleotide_diverstiy': population_nucleotide_diversity(df_snv, ref_seq_length)})

    # mean position-wise Shannon entropy
    out_dict.update({'mean_position_shannon': mean_pos_Shannon_entropy(df_snv,ref_seq_length)})
    df_diveristy = pd.DataFrame(list(out_dict.items()))
    df_diveristy.to_csv(output_dir + 'diversity_measures.csv', index=False, header=False)

    # position-wise Shannon entropy
    pos_shannon_dict={}
    for i in range(ref_seq_length):
        if position_Shannon_entropy(df_snv,i)!= 0:
            pos_shannon_dict.update({i: position_Shannon_entropy(df_snv,i)})
    df_pos_shannon = pd.DataFrame(list(pos_shannon_dict.items()), columns=['position', 'Shannon_entropy'])
    df_pos_shannon.to_csv(output_dir + 'position_shannon_entropy.csv', index=False, header=False)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
