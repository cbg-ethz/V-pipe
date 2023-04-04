#!/usr/bin/env python3
import argparse
import pandas as pd
import regex
import ruamel.yaml
import datetime
import numpy as np
import os
import sys

regex.DEFAULT_VERSION = regex.VERSION1

""" This script is the parse_samname function from mut_table notebook"""


def parse_args():
    """Set up of command-line arguments parsing
    Arguments:
        sampmle_file: input tsv file of samples
        -l: wwtp codes
        -c: config file including regex and custom configs
        -o: output file
        -n/-0/--lower: help formatting corner cases
    """

    parser = argparse.ArgumentParser(description="Process information of sample data")

    parser.add_argument(
        "sample_file",
        type=str,
        help="TSV file of samples",
    )

    parser.add_argument(
        "-l",
        "--locations",
        help="maps locations (wastewater treatment plant) short code to full name",
        default=None,
        dest="loc_file",
        required=False,
    )

    parser.add_argument(
        "-c",
        "--regex-config",
        help='yaml file for regex, MUST provide groups "location", "year", "month", and "day"',
        default=None,
        dest="config_file",
        required=False,
    )

    parser.add_argument(
        "-n",
        "--no-fallback",
        help="do not fallback to copy if regex do not match",
        action="store_false",
        dest="fallback",
    )

    parser.add_argument(
        "-0",
        "--no-strip0",
        help='do not strip leading zeroes from numbers (i.e.: treat "01" and "1" as distinct)',
        action="store_false",
        dest="strip0",
    )

    parser.add_argument(
        "--lower",
        help='ignore case (i.e.: treat "cOdE" is if it were "code")',
        action="store_false",
        dest="lower",
    )

    parser.add_argument(
        "-o", "--output", help="output file", dest="output", required=True
    )

    parser.add_argument(
        "--out-locations",
        help="write list of found locations",
        dest="out_loc",
        default=None,
    )

    return parser.parse_args()


def parse_samname(
    sample,
    batch,
    rxsam=None,
    rxbat=None,
    datefmt=None,
    fallback=True,
    strip0=True,
    lower=True,
):
    """Parse tsv file of samples into: date, location_code"""

    # use the regex (if provided) to interpret sample and batch
    m_sam = {}
    m_bat = {}

    if rxsam:
        match_sam = rxsam.search(sample)
        if match_sam:
            m_sam = match_sam.groupdict()
    if rxbat:
        match_batch = rxbat.search(batch)
        if match_batch:
            m_bat = match_batch.groupdict()

    # search order:
    # - matched in sample name?
    # - matched in batch ?
    # - fall back to simple copy

    # location
    location_code = (
        m_sam.get("location", None)
        or m_bat.get("location", None)
        or (sample if fallback else "")
    )
    location_extra = m_sam.get("location_extra", None) or m_bat.get(
        "location_extra", None
    )
    if location_extra:
        location_code += "_" + location_extra

    # date
    date = (
        m_sam.get("date", None)
        or m_bat.get("date", None)
        or (batch if fallback else "")
    )

    # prefered date: explicit year/month/day
    y = m_sam.get("year", None) or m_bat.get("year", None)
    if y is not None:
        y = int(y)
        if y < 100:
            y += 2000
    m = m_sam.get("month", None) or m_bat.get("month", None)
    if m is not None:
        m = int(m)
    d = m_sam.get("day", None) or m_bat.get("day", None)
    if d is not None:
        d = int(d)

    if y and m and d:
        date = datetime.datetime(y, m, d).strftime("%Y-%m-%d")
    # fallback date:
    elif datefmt and date:
        try:
            # try as-is
            date = datetime.datetime.strptime(date, datefmt).date().strftime("%Y-%m-%d")
        except ValueError as v:
            # drop any extra characters
            if len(v.args) > 0 and v.args[0].startswith("unconverted data remains: "):
                print(f"Warning: {v.args[0]}", file=sys.stderr)
                date = date[: -(len(v.args[0]) - 26)]
                date = (
                    datetime.datetime.strptime(date, datefmt)
                    .date()
                    .strftime("%Y-%m-%d")
                )
            else:
                raise

    # normalization/sanitation:
    # - remove leading zeros, only from numbers
    if strip0 and location_code.isdigit():
        location_code = location_code.lstrip("0") or "0"
    # - lowercase
    if lower:
        location_code = location_code.lower()

    return location_code, date


def write_tsv(sample_name, sample_info, out_folder):
    """Write tsv file of sample information"""
    info_file = f"{out_folder}/{sample_name}.tsv"
    if not os.path.exists(info_file):
        with open(info_file, "w") as f:
            f.write("%s %s %s %s \n" % sample_info)
    else:
        with open(info_file, "a") as f:
            f.write("%s %s %s %s \n" % sample_info)


def main():
    args = parse_args()
    assert os.path.exists(
        args.sample_file
    ), f"Cannot find input file {args.sample_file}"

    # load regexes
    rxsam = None
    rxbat = None
    datefmt = None
    if args.config_file:
        with open(args.config_file, "r") as stream:
            try:
                reg_data = ruamel.yaml.load(stream, Loader=ruamel.yaml.Loader)
            except yaml.YAMLError as exc:
                print(exc)
                sys.exit(1)
        if "sample" in reg_data:
            rxsam = regex.compile(reg_data.get("sample"))
        if "batch" in reg_data:
            rxbat = regex.compile(reg_data.get("batch"))
        datefmt = reg_data.get("datefmt")

    # load location full names
    locations = None
    if args.loc_file:
        locations = pd.read_csv(args.loc_file, sep="\t")
        if args.strip0:
            locations["code"] = locations["code"].apply(
                lambda s: (s.lstrip("0") or "0") if s.isdigit() else s
            )
        if args.lower:
            locations["code"] = locations["code"].str.lower()

        locations = locations.set_index("code")

    # samples tsv input
    samples_info = pd.read_csv(
        args.sample_file, sep="\t", names=["sample", "batch", "reads", "proto"]
    )

    # add location code and date column
    if rxsam or rxbat or datefmt:
        # extract information using regexs
        samples_info[["location_code", "date"]] = samples_info[
            ["sample", "batch"]
        ].apply(
            lambda r: parse_samname(
                sample=r["sample"],
                batch=r.batch,
                rxsam=rxsam,
                rxbat=rxbat,
                datefmt=datefmt,
                fallback=args.fallback,
                strip0=args.strip0,
                lower=args.lower,
            ),
            axis=1,
            result_type="expand",
        )
    else:
        # simple copy: 'sample' name is 'location' and 'batch' is the 'date'
        samples_info[["location_code", "date"]] = samples_info[["sample", "batch"]]

    # add location full name
    if locations is not None:
        samples_info = samples_info.merge(
            locations, how="left", left_on="location_code", right_index=True
        )

    samples_info.set_index(["sample", "batch"]).to_csv(
        args.output, sep="\t", compression={"method": "infer"}
    )

    if args.out_loc:
        with open(args.out_loc, "w") as outf:
            ruamel.yaml.round_trip_dump(
                {
                    "locations_list": list(
                        set(
                            samples_info[
                                "location" if locations is not None else "location_code"
                            ].unique()
                        )
                        - {"", np.nan}
                    )
                },
                outf,
            )


if __name__ == "__main__":
    main()
