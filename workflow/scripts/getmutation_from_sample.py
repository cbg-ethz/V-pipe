#!/usr/bin/env python3
"""
Script retrieving the coverage of mutations from a single sample
"""
import pandas as pd
import numpy as np
import os
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

    requiredNamed = parser.add_argument_group("arguments used for V-pipe integration")
    requiredNamed.add_argument(
        "-s",
        "--samplename",
        required=False,
        type=str,
        dest="sample",
        help="sample_name as found in the V-pipe working directory",
    )
    requiredNamed.add_argument(
        "-b",
        "--batch",
        required=False,
        type=str,
        dest="batch",
        help="Name of the batch the sample is part of",
    )
    parser.add_argument(
        "-o",
        "--output",
        "--outname",
        required=False,
        type=str,
        dest="outname",
        help="Filename of the final output table. If not provided, it defaults to <samplename>_mutations.txt",
    )
    parser.add_argument(
        "-m",
        "--mutationtable",
        required=False,
        default="mutlist.txt",
        type=str,
        dest="muttable",
        help="Mutations helper table",
    )
    parser.add_argument(
        "-a",
        "--based",
        required=False,
        default=1,
        type=int,
        dest="base",
        help="Are the positions in the tsv 0-based or 1-based?",
    )
    parser.add_argument(
        "basecnt",
        metavar="BASECOUNT",
        type=str,
        help="basecount TSV table produced by V-pipe to search for mutations",
    )

    return parser.parse_args()


#####


def getmutation_from_sample(basecnt, tsvbase, mut):
    # warning that table is *tsvbase*-based
    basecount = (
        pd.read_csv(
            basecnt,
            sep="\t",
            header=[0, 1],
            index_col=[0, 1],
        )
        .droplevel("ref")
        .T.droplevel("sample")
        .T
    )
    # total coverage
    basecount["cov"] = basecount.apply(sum, axis=1)
    # look mutations per position
    return pd.DataFrame(
        data=mut.apply(
            lambda x: pd.concat(
                [
                    pd.Series(
                        [
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
    )


###
"""
Helper functions
"""


def build_outname(args):
    if args.outname is None:
        args.outname = f"{args.sample}_{args.batch}_mutations.txt"
    return args


###


def main():
    args = parse_args()
    args = build_outname(args)

    # list of mutations to search
    mut = pd.read_csv(args.muttable, sep="\t").astype({"position": "int"})

    # seach them!
    assert os.path.exists(args.basecnt), f"Cannot find input file {args.basecnt}"
    table = getmutation_from_sample(basecnt=args.basecnt, tsvbase=args.base, mut=mut)
    assert table.shape[0] > 0, "Generated an empty mutation table!"

    idx = []
    # add extra columns
    if args.sample:
        table["sample"] = args.sample
        idx += ["sample"]
    if args.batch:
        table["batch"] = args.batch
        idx += ["batch"]

    # set index
    idx += ["pos"]
    table.set_index(idx, inplace=True)

    print(args.outname)
    table.to_csv(args.outname, sep="\t", compression={"method": "infer"})


if __name__ == "__main__":
    main()
