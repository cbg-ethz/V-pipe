from pathlib import Path
import pandas as pd

def runtime(benchmark_list, fname_out):
    # gather benchmark information
    tmp = []
    for fname in benchmark_list:
        parts = str(fname).split("/")
        if len(parts) == 7:
            _, _, params, method, _, replicate, _ = parts
        elif len(parts) == 8:  # for multi workflow
            _, _, _, params, method, _, replicate, _ = parts

        df_tmp = pd.read_csv(fname, sep="\t")
        df_tmp["method"] = method
        df_tmp["params"] = params
        df_tmp["replicate"] = replicate

        tmp.append(df_tmp)
    df_bench = pd.concat(tmp).replace("-", pd.NA)

    # plot
    df_long = (
        pd.melt(df_bench, id_vars=["method", "params", "replicate"])
        .assign(params=lambda x: x["params"].str.replace("_", "\n"))
        .query("variable == 's'")
    )
    df_long.to_csv(fname_out)


def main(benchmark_list, fname_out):

    runtime(benchmark_list, fname_out)


if __name__ == "__main__":
    main(
        [Path(e) for e in snakemake.input.benchmark_list],
        snakemake.output.fname_out,
    )
