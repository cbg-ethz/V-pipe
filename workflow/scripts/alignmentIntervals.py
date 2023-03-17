#!/usr/bin/env python3

import os
import argparse
from math import floor
from shutil import copyfile

DBG = True if os.environ.get("DBG") is not None else False


def parse_args():
    """Set up the parsing of command-line arguments"""

    parser = argparse.ArgumentParser(
        description="Benchmark: compute the union/intersection among coverage "
        "intervals to compare performance between aligner or SNV "
        "caller",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    requiredNamed = parser.add_argument_group("required named arguments")
    requiredNamed.add_argument(
        "-d",
        required=True,
        default=None,
        metavar="DIR1,DIR2,...",
        dest="directories",
        help="Sub-directories where results from vpipeBench are stored",
    )
    parser.add_argument(
        "-r",
        required=False,
        default=None,
        metavar="FASTA",
        dest="reference_file",
        type=str,
        help="Referene sequence",
    )
    parser.add_argument(
        "-wl",
        required=False,
        default=None,
        metavar="len1,len2,...",
        dest="window_len",
        type=str,
        help="Window length used by ShoRAH to construct overlapping windows",
    )
    parser.add_argument(
        "-ws",
        required=False,
        default=None,
        metavar="shift1,shift2,...",
        dest="window_shift",
        type=str,
        help="Length of the window shift used by ShoRAH to construct "
        "overlapping windows",
    )
    parser.add_argument(
        "--caller",
        required=False,
        default=None,
        action="store_true",
        dest="comparison_callers",
        help="Indicate whether to construct alignment intervals for comparing "
        "mutation callers",
    )
    parser.add_argument(
        "-o",
        required=False,
        default=None,
        metavar="DIR",
        dest="outdir",
        type=str,
        help="Output directory",
    )

    return parser.parse_args()


def read_fasta(fasta_file):
    "Parse fasta file"
    count = 0
    headers = []
    sequences = []
    aux = []
    with open(fasta_file, "r") as infile:
        for line in infile:
            record = line.rstrip()
            if record and record[0] == ">":
                headers.append(record[1:])
                if count > 0:
                    sequences.append("".join(aux))
                    aux = []
            else:
                aux.append(record)
            count += 1

    sequences.append("".join(aux))
    return headers, sequences


def parse_region(region, shorah, reference_file, win_len, win_shift):
    def parse_region_aux(region, shorah, reference_file, win_len, win_shift):
        aux = region.split(":")
        reference_name = aux[0]
        aux = aux[1].split("-")
        start = int(aux[0])
        end = int(aux[1])

        if shorah:
            header, reference = read_fasta(reference_file)
            reference_len = len(reference[0])
            # ShoRAH was used for SNV calling
            # Assuming 3 windows were used for SNV calling, identify region
            # that is spanned by at least 2 windows (below, using
            # 0-based indexing and closed intervals)
            start_ = max(0, start - win_len - 1)
            end_ = min(reference_len, end + win_len)
            num_windows = floor((end_ - (start_ + win_len - 1)) / win_shift) + 1
            offset = 2 * win_shift

            start = max(0, start - offset - 1)
            # In order to identify the region which is covered by at least
            # two windows, add to the end of the first window the
            # increment multiply by the number of windows - 2 (i.e.,
            # discarding last window). In this case assuming half-open
            # interval [start, end)
            end = min(reference_len, start_ + win_len + (num_windows - 2) * win_shift)

        return reference_name, start, end

    aux = region.split(",")
    reference_name = []
    interval = []
    for s in aux:
        ref_, start_, end_ = parse_region_aux(
            s, shorah, reference_file, win_len, win_shift
        )
        reference_name.append(ref_)
        interval.append((start_, end_))

    return reference_name, interval


def parse_intervals(
    input_file,
    union_intervals,
    intersection_intervals,
    shorah,
    reference_file,
    window_len,
    window_shift,
):
    def merge_regions(intervals_new, intervals_old):
        """
        Merge intervals by:
        (1) sorting the intervals in decreasing order, and
        (2) iterate through intervals, comparing the current interval with the
            last one. If the interval overlap, merge both intervals into one
        """
        intervals = intervals_old + intervals_new
        intervals_sorted = sorted(intervals, key=lambda tup: tup[0])
        intervals_merged = []

        for interval in intervals_sorted:
            if not intervals_merged:
                intervals_merged.append(interval)
            else:
                last = intervals_merged[-1]
                # Test for overlap between current interval and the last one.
                # We know that last[0] <= interval[0]
                if interval[0] <= last[1]:
                    upper_bound = max(last[1], interval[1])
                    intervals_merged[-1] = (last[0], upper_bound)
                else:
                    intervals_merged.append(interval)

        return intervals_merged

    def find_intersect(intervals_new, intervals_old):
        intervals_intersected = []

        for i_new in intervals_new:
            for i_old in intervals_old:
                lower_bound = max(i_old[0], i_new[0])
                upper_bound = min(i_old[1], i_new[1])

                if lower_bound < upper_bound:
                    intervals_intersected.append((lower_bound, upper_bound))

        return intervals_intersected

    idx = 0
    with open(input_file, "r") as infile:
        for line in infile:
            record = line.rstrip().split("\t")
            key = record[0]
            if len(record) > 1:
                win_len = int(window_len[idx]) if window_len else None
                win_shift = int(window_shift[idx]) if window_shift else None
                reference_name, intervals = parse_region(
                    record[1], shorah, reference_file, win_len, win_shift
                )
            else:
                reference_name, intervals = None, None
            idx += 1

            if key in union_intervals:
                # Update region if it was empty or if has to be enlarged
                if not union_intervals[key][1]:
                    if DBG:
                        print(f"Assigning {key}: {reference_name}, {intervals}")
                    union_intervals[key] = [reference_name, intervals]
                elif intervals:
                    intervals_new = merge_regions(intervals, union_intervals[key][1])
                    # It assumes a unique reference sequence
                    union_intervals[key][0] = [
                        reference_name[0] for i in range(len(intervals_new))
                    ]
                    union_intervals[key][1] = intervals_new
                    if DBG:
                        print(
                            f"Merging intervals for {key}: "
                            f"{union_intervals[key][0]}, "
                            f"{union_intervals[key][1]}"
                        )
            else:
                if DBG:
                    print(f"Assigning {key}: {reference_name}, {intervals}")
                union_intervals[key] = [reference_name, intervals]

            if key in intersection_intervals:
                if intersection_intervals[key][1] and intervals:
                    intervals_new = find_intersect(
                        intervals, intersection_intervals[key][1]
                    )
                    if intervals_new:
                        # It assumes a unique reference sequence
                        intersection_intervals[key][0] = [
                            reference_name[0] for i in range(len(intervals_new))
                        ]
                    intersection_intervals[key][1] = intervals_new
                elif not intervals:
                    # One empty interval is enough to set the intersect to the
                    # empty set
                    intersection_intervals[key] = None, None
            else:
                intersection_intervals[key] = [reference_name, intervals]

    return union_intervals, intersection_intervals


def read_coverage_file(input_file, reference_file, window_len, window_shift):
    out_dict = {}
    idx = 0
    with open(input_file, "r") as infile:
        for line in infile:
            record = line.rstrip().split("\t")
            key = record[0]
            if len(record) > 1:
                win_len = int(window_len[idx]) if window_len else None
                win_shift = int(window_shift[idx]) if window_shift else None
                reference_name, intervals = parse_region(
                    record[1], True, reference_file, win_len, win_shift
                )
            idx += 1
            out_dict[key] = [reference_name, intervals]

    return out_dict


def write_output(output_file, intervals_dir):
    with open(output_file, "w") as outfile:
        for key, value in intervals_dir.items():
            if value[0] and value[1]:
                str_ = []
                for idx, ref in enumerate(value[0]):
                    region = value[1][idx]
                    str_.append(f"{ref}:{region[0]}-{region[1]}")
                str_ = ",".join(str_)
                outfile.write(f"{key}\t{str_}\n")
            else:
                outfile.write(f"{key}\t\n")


def main():
    args = parse_args()

    union_dir_shorah = {}
    intersection_dir_shorah = {}

    union_dir = {}
    intersection_dir = {}

    if args.window_len:
        win_len = args.window_len.split(",")
    if args.window_shift:
        win_shift = args.window_shift.split(",")

    outdir = args.outdir if args.outdir is not None else os.getcwd()

    # Aggregarting results
    directories = args.directories.split(",")
    for indir in directories:
        if args.comparison_callers:
            # For the comparison of tools for mutation calling, the alignment
            # remains a constant. However, ShoRAH only reports SNVs that are
            # covered by at least two of the overlapping windows.
            # "union", in this context, results in using all positions above
            # certain coverage
            # "intersect", accounts for the overlapping windows.
            aligner = indir.split("-")[0]
            if "shorah" in indir:
                input_file = os.path.join(indir, "variants", "coverage_intervals.tsv")
                intersect_caller = read_coverage_file(
                    input_file, args.reference_file, win_len, win_shift
                )
                output_file = os.path.join(
                    outdir, f"intersect_coverage_intervals_{aligner}.tsv"
                )
                write_output(output_file, intersect_caller)
            else:
                src = os.path.join(
                    os.getcwd(), indir, "stats", "coverage_intervals.tsv"
                )
                output_file = os.path.join(
                    outdir, f"union_coverage_intervals_{aligner}.tsv"
                )
                copyfile(src, output_file)
        else:
            if "shorah" in indir:
                input_file = os.path.join(indir, "variants", "coverage_intervals.tsv")
                union_dir_shorah, intersection_dir_shorah = parse_intervals(
                    input_file,
                    union_dir_shorah,
                    intersection_dir_shorah,
                    True,
                    args.reference_file,
                    win_len,
                    win_shift,
                )
            else:
                input_file = os.path.join(indir, "stats", "coverage_intervals.tsv")
                union_dir, intersection_dir = parse_intervals(
                    input_file,
                    union_dir,
                    intersection_dir,
                    False,
                    None,
                    None,
                    None,
                )

    if not args.comparison_callers:
        output_file = os.path.join(outdir, "union_coverage_intervals_ShoRAH.tsv")
        write_output(output_file, union_dir_shorah)

        output_file = os.path.join(outdir, "intersect_coverage_intervals_ShoRAH.tsv")
        write_output(output_file, intersection_dir_shorah)

        output_file = os.path.join(outdir, "union_coverage_intervals.tsv")
        write_output(output_file, union_dir)

        output_file = os.path.join(outdir, "intersect_coverage_intervals.tsv")
        write_output(output_file, intersection_dir)


if __name__ == "__main__":
    main()
