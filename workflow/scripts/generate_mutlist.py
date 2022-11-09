#!/usr/bin/env python3
import pandas as pd
import strictyaml
import re
import glob
import os
import sys
import click

rxmutdec = re.compile(
    r"^(?:(?:(?:(?P<ref>[ATCG]+)\>)?(?P<mut>[ATCG]+))|(?P<del>[\-]+)|(?:[\+](?P<ins>[ATGC]+)))$"
)

schema = {
    "position": [],
    "reference": [],
    "variant": [],
}

types = {"position": "int"}


def load_voc_yaml(yp):
    with open(yp, "r") as yf:
        yam = strictyaml.dirty_load(yf.read(), allow_flow_style=True).data
    muts = pd.DataFrame(
        data=schema
        | {
            yam["variant"]["short"]: [],
        }
    ).astype(types)
    for c in ["mut", "revert", "extra", "shared", "subset"]:
        # all categories (we don't care, we will compare accross samples)
        if c not in yam:
            continue

        for pos, mutstr in yam[c].items():
            if not (res := rxmutdec.match(mutstr)):
                print(f"{yp}:{pos} cannot parse {mutstr}")
                continue
            match = res.groupdict()
            if match["ins"]:
                print(f"{yp}:{pos} insertions not supported (yet): {match['ins']}")
                continue
            elif match["mut"]:
                for i in range(len(match["mut"])):
                    muts = pd.concat(
                        [
                            muts,
                            pd.DataFrame.from_records(
                                [
                                    {
                                        "position": int(pos) + i,
                                        "reference": (
                                            match["ref"][i]
                                            if match["ref"] and i < len(match["ref"])
                                            else "N"
                                        ),
                                        "variant": match["mut"][i],
                                        yam["variant"]["short"]: c,
                                    }
                                ]
                            ),
                        ]
                    )
            elif match["del"]:
                # TODO this is wrong and will be fixed in ShoRAH
                for i in range(len(match["del"])):
                    muts = pd.concat(
                        [
                            muts,
                            pd.DataFrame.from_records(
                                [
                                    {
                                        "position": int(pos) + i,
                                        "reference": (
                                            match["ref"][i]
                                            if match["ref"] and i < len(match["ref"])
                                            else "N"
                                        ),
                                        "variant": "-",
                                        yam["variant"]["short"]: c,
                                    }
                                ]
                            ),
                        ]
                    )
    return muts


@click.command(
    help="Generate the mutlist used when looking for variant using variant signatures",
    # epilog="",
)
@click.option(
    "--output",
    "--out",
    "-o",
    metavar="TSV",
    required=False,
    default="mutlist.tsv",
    type=click.Path(),
    help="Write results to this output TSV instead of 'mutlist.tsv'",
)
@click.option(
    "--genes",
    "-g",
    metavar="GFF",
    required=False,
    default=None,
    type=click.Path(exists=True),
    help="Add 'gene' column to table",
)
@click.option(
    "--voc-dir",
    "-d",
    metavar="PATH",
    multiple=True,
    required=False,
    default=[],
    type=str,
    help="scan directory for additional voc YAML files",
)
@click.option(
    "--verbose/--no-verbose",
    "-v/-V",
    required=False,
    default=False,
    type=bool,
    help="verbose (dumps table on terminal)",
)
@click.argument("vocs", metavar="VOC_YAML", nargs=-1, type=click.Path(exists=True))
def generate_mutlsit(output, genes, voc_dir, vocs, verbose):
    # scan directories for addition voc YAMLS
    if len(voc_dir):
        vocs = list(vocs) + [
            yp for vd in voc_dir for yp in glob.glob(os.path.join(vd, "*.yaml"))
        ]  # auto skips .hidden

    assert len(
        vocs
    ), f"at least provide some voc YAML files, by listing them directly or either by scanning directories with option '--voc-dir'."

    vartable = pd.DataFrame(data=schema).astype({"position": "int"})
    for yp in set(vocs):
        print(yp)
        assert os.path.exists(yp), f"cannot find {yp}"

        muts = load_voc_yaml(yp)

        vartable = vartable.merge(
            how="outer", right=muts, copy=False, sort=True
        )  # .fillna('')

    if genes:
        from BCBio import GFF

        # place an empty column right after the standard columns
        if not "gene" in vartable.columns:
            vartable.insert(len(schema), "gene", [""] * len(vartable.index))
        with open(genes) as gf:
            for record in GFF.parse(gf):
                for feature in record.features:
                    if feature.type == "gene":
                        mask = (int(feature.location.end) >= vartable["position"]) & (
                            vartable["position"] >= int(feature.location.start)
                        )
                        vartable.loc[mask, "gene"] = feature.qualifiers.get(
                            "Name", [feature.id]
                        )[0]

    if verbose:
        with pd.option_context(
            "display.max_rows", None
        ):  # , 'display.max_columns', None):
            print(vartable, file=sys.stderr)  # .sort_values('position'))

    vartable.to_csv(
        output, sep="\t", compression={"method": "infer"}, index=False, na_rep="NA"
    )


if __name__ == "__main__":
    generate_mutlsit()
