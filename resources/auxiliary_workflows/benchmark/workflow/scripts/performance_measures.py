from pathlib import Path


def main(vcf_list, groundtruth_list, dname_out):
    # TODO: do something
    dname_out.mkdir(parents=True)

    for fname_vcf, fname_groundtruth in zip(vcf_list, groundtruth_list):
        print(fname_vcf, fname_groundtruth)


if __name__ == "__main__":
    main(
        [Path(e) for e in snakemake.input.vcf_list],
        [Path(e) for e in snakemake.input.groundtruth_list],
        Path(snakemake.output.dname_out),
    )
