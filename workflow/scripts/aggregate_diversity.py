#!/usr/bin/env python3
"""
Script aggregating diversity measures from all samples into one file.
"""
import pandas as pd


def main(in_fnames_diversity, in_fnames_shannon, out_diversity_csv, out_shannon_csv):
    merged_div_csv = pd.concat(
        [pd.read_csv(path_div) for path_div in in_fnames_diversity]
    )
    merged_div_csv.to_csv(out_diversity_csv)

    merged_shan_csv = pd.concat(
        [pd.read_csv(path_shan) for path_shan in in_fnames_shannon]
    )
    merged_shan_csv.to_csv(out_shannon_csv)


if __name__ == "__main__":
    main(
        snakemake.input.fnames_diversity,
        snakemake.input.fnames_shannon,
        snakemake.output.diversity_csv,
        snakemake.output.shannon_csv,
    )
