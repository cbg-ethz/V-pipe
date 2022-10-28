#!/usr/bin/env python3
"""
Script retrieving the coverage of mutations from a single sample
"""
import pandas as pd
import numpy as np
import datetime
import os
import re
import argparse
import sys

__author__ = "Matteo Carrara"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"

def parse_args():
    """Set up the parsing of command-line arguments"""

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    requiredNamed = parser.add_argument_group("required named arguments")
    requiredNamed.add_argument(
        '-s',
        '--samplename',
        required=True,
        type=str,
        dest="sample",
        help="sample_name as found in the V-pipe working directory"
    )
    requiredNamed.add_argument(
        '-b',
        '--batch',
        required=True,
        type=str,
        dest="batch",
        help="Name of the batch the sample is part of"
    )
    requiredNamed.add_argument(
        '-p',
        '--proto',
        required=True,
        type=str,
        dest="proto",
        help="Amplification protocol used. Format is <v#>"
    )
    requiredNamed.add_argument(
        '-d',
        '--samplesdir',
        required=True,
        type=str,
        dest="samples_dir",
        help="Samples directory"
    )
    parser.add_argument(
        '-o',
        '--outname',
        required=False,
        type=str,
        dest="outname",
        help="Filename of the final output table. If not provided, it defaults to <samplename>_mutations.txt"
    )
    parser.add_argument(
        '-m',
        '--mutationtable',
        required=False,
        default="mutlist.txt",
        type=str,
        dest="muttable",
        help="Mutations helper table"
    )
    parser.add_argument(
        '-n',
        '--parsedname',
        required=False,
        default="parsed_sample_name.tsv",
        type=str,
        dest="parsed_name_table",
        help="3-field table resulting from the parsing of the sample name"
    )
    parser.add_argument(
        '-a',
        '--based',
        required=False,
        default=1,
        type=int,
        dest="base",
        help="Are the positions in the tsv 0-based or 1-based?"
    )
    
    return parser.parse_args()
#####

def getmutation_from_sample(tsam, tbat, tproto, mut, samples_dir, date, plantcode, plantname, tsvbase):
    # warning that table is *tsvbase*-based
    basecount = (
        pd.read_csv(
            f"{samples_dir}/{tsam}/{tbat}/alignments/basecnt.tsv.gz",
            sep="\t",
            header=[0, 1],
            index_col=[0, 1],
        )
        .droplevel("ref")
        .T.droplevel("sample")
        .T
    )
    basecount["cov"] = basecount.apply(sum, axis=1)
    r = pd.DataFrame(
        data=mut.apply(
            lambda x: pd.concat(
                [
                    pd.Series(
                        [
                            tsam,
                            tbat,
                            tproto,
                            date,
                            plantcode,
                            plantname,
                            x.gene,
                            x.position,
                            x.variant,
                            # -1 : 1-based to 0-based
                            basecount.loc[x.position - (1 - tsvbase)]["cov"],
                            basecount.loc[x.position - (1 - tsvbase)][x.variant],
                            basecount.loc[x.position - (1 - tsvbase)][x.variant]
                            / basecount.loc[x.position - (1 - tsvbase)]["cov"]
                            if basecount.loc[x.position - (1 - tsvbase)]["cov"]
                            else np.nan,
                        ],
                        index=[
                            "sample",
                            "batch",
                            "proto",
                            "date",
                            "plantcode",
                            "plantname",
                            "gene",
                            "pos",
                            "base",
                            "cov",
                            "var",
                            "frac",
                        ],
                    ),
                    pd.Series(x[4:]),
                ]
            ),
            axis=1,
        )
    ).set_index(["sample", "batch", "pos"])
    # testing
    #     if b:
    #         print(r)
    return r

###
"""
Helper functions
"""

def build_outname(args):
    if args.outname is None:
        args.outname = f"{args.sample}_{args.batch}_mutations.txt"
    return args

def load_parser_table(args):
    parsed_name = pd.read_csv(args.parsed_name_table, sep="\t")
    parsed_name = parsed_name[parsed_name["batch"]==args.batch]
    return parsed_name

###

def main():
    args = parse_args()
    args = build_outname(args)
    mut = pd.read_csv(args.muttable, sep="\t").astype({"position": "int"})
    #parsed_name = pd.read_csv(args.parsed_name_table, sep="\t")
    parsed_name = load_parser_table(args)
    if parsed_name.shape[0] == 0:
        print("Provided batch does not match any batch in the name parser table!")
        sys.exit()
    date = parsed_name["date"].values[0]
    plantcode = parsed_name["plantcode"].values[0]
    plantname = parsed_name["plantname"].values[0]
    batch = parsed_name["batch"].values[0]
    table = getmutation_from_sample(args.sample, batch, args.proto, mut, args.samples_dir, date, plantcode, plantname, args.base)
    if table.shape[0] == 0:
        print('Generated an empty mutation table!')
        sys.exit()
    print(args.outname)
    table.to_csv(args.outname, sep="\t") 

if __name__ == "__main__":
    main()
