#!/usr/bin/env python3
"""
Script retrieving all single-sample mutation lists and merging them together in a single tally table
"""

import pandas as pd
import numpy as np
import datetime
import os
import argparse

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
        '--sampleslist',
        required=True,
        type=str,
        dest="sampleslist",
        help="table containing the list of samples to merge together"
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
        '-i',
        '--inputname',
        required=False,
        default="_mutations.txt",
        type=str,
        dest="inputname",
        help="Name structure of the input name to use. If present, used to create the structure <sample_name><inputname>. If absent, assumed to be <sample_name>_mutations.txt"
    )
    parser.add_argument(
        '-f',
        '--failedsamplesfile',
        required=False,
        default="failed_samples.txt",
        type=str,
        dest="failed_samples_filename",
        help="Filename (and, optionally, full path) where to store the list of provided samples for which a mutation table could not be found"
    )
    parser.add_argument(
        '-o',
        '--outname',
        required=False,
        default="tallymut.tsv",
        type=str,
        dest="outname",
        help="Filename (and, optionally, full path) where to store the full table of mutations for all samples"
    )

    return parser.parse_args()

### 
"""
Helper functions
"""

def build_input_path(basepath, sname, batch, inputname, failed_samples_filename):
    input_filename = f"{sname}_{batch}{inputname}"
    input_path = f"{basepath}/{sname}/{batch}/{input_filename}"
    try:
        table = pd.read_csv(input_path, sep="\t")
    except FileNotFoundError:
        print(f"Mutation file unavailable for listed sample {sname} in batch {batch}. Please check if directory {input_path} is correct")
        with open(failed_samples_filename, 'a') as f:
            current = datetime.datetime.now()
            f.write(f"{current} {input_path}\n")
            return None
    return table

###

def main():
    args = parse_args()
    sample_tbl = pd.read_csv(args.sampleslist, sep="\t", header=None)
    all_input_files = [build_input_path(args.samples_dir, sname, batch, args.inputname, args.failed_samples_filename) for sname,batch in zip(sample_tbl[0], sample_tbl[1])]
    full_table = pd.concat(all_input_files)
    full_table.to_csv(args.outname, sep="\t", index=False)

if __name__ == "__main__":
    main()

