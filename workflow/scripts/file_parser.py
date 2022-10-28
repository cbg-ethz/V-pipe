import argparse
import pandas as pd
import regex
import yaml
import datetime
import numpy as np
import os

regex.DEFAULT_VERSION = regex.VERSION1

""" This script is the parse_samname function from mut_table notebook"""


def parse_args():
    """ Set up of command-line arguments parsing
    Arguments:
        -i: input tsv file of samples
        -p: wwtp codes
        -c: config file including regex and custom configs
        -o: output folder for separate tsv files
    """

    parser = argparse.ArgumentParser(description="Process information of sample data")

    parser.add_argument(
        '-i',
        help='TSV file of samples',
        dest='sample_file',
        required=True
    )

    parser.add_argument(
        '-p',
        help='WW plants code',
        dest='plants_file',
        required=True
    )

    parser.add_argument(
        '-c',
        help='yaml file for regex',
        dest='config_file',
        required=True
    )

    parser.add_argument(
        '-o',
        help='output folder',
        dest='out_fold',
        required=True
    )

    return parser.parse_args()


def parse_samname(sample_info, reg_info, ww_plants, name_conv):
    """Parse tsv file of samples into: date, plantcode, plantname, batch"""

    date = plantcode = plantname = np.nan
    tsam = sample_info['sample']
    tbam = sample_info['batch']

    if name_conv == 'default':
        rxname = regex.compile(reg_info['regex']['regex_default'])
        match = rxname.search(tsam)
        if match:
            m = match.groupdict()
            # print(m)
            if m['default']:
                if m["month"] and m["day"]:
                    date = datetime.datetime(
                        int(m["year"]), int(m["month"]), int(m["day"])
                    ).strftime("%Y-%m-%d")
                plantcode = int(m["plant"])
                plantname = (
                    ww_plants.at[plantcode, "Plant"] if plantcode in ww_plants.index else ""
                )
            elif m["other"]:
                if m["KLZH"]:
                    # print('>>>>>>>>>>', tsam, m)
                    date = (
                        datetime.datetime.strptime(m["date"], "%y%m%d").date().strftime("%Y-%m-%d")
                    )
                    if not m["KLZHsuffix"]:  # avoid _Promega and _2
                        plantname = "Kanton Zürich"
                        plantcode = 90
                    else:
                        plantname = "Kanton Zürich/Promega"
                        plantcode = 91
                elif m["BA"]:
                    if tsam in reg_info['basel_patchmap'].keys():
                        date = reg_info['basel_patchmap'][tsam]
                    elif m["date"]:
                        date = (
                            datetime.datetime.strptime(m["date"], "%Y-%m-%d").date().strftime("%Y-%m-%d")
                        )
                    plantname = "Basel (catchment area ARA Basel)"
                    plantcode = 92

    return date, plantcode, plantname, tbam


def write_tsv(sample_name, sample_info, out_folder):
    """Write tsv file of sample information"""
    info_file = f'{out_folder}/{sample_name}.tsv'
    if not os.path.exists(info_file):
        with open(info_file, 'w') as f:
            f.write('%s %s %s %s \n' % sample_info)
    else:
        with open(info_file, 'a') as f:
            f.write('%s %s %s %s \n' % sample_info)


def main():
    args = parse_args()
    if not [x for x in (args.plants_file, args.sample_file, args.config_file) if x is None]:

        with open(args.config_file, "r") as stream:
            try:
                reg_data = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)

        plants = pd.read_csv(args.plants_file, sep='\t')
        plants = plants.set_index('Code')
        samples_info = pd.read_csv(args.sample_file, sep='\t', names=["sample", "batch", "reads", "proto"])
        c = 0
        for i, samp in samples_info.iterrows():
            info = parse_samname(samp, reg_data, plants, 'default')
            write_tsv(samp['sample'], info, args.out_fold)


if __name__ == '__main__':
    main()
