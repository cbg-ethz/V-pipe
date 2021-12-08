from pathlib import Path


def main(vcf_list, dname_out):
    # TODO: do something
    dname_out.mkdir(parents=True)


if __name__ == "__main__":
    main([Path(e) for e in snakemake.input.vcf_list], Path(snakemake.output.dname_out))
