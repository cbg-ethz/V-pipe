"""Amend genotype field to contain all variants."""

import sys

from cyvcf2 import VCF, Writer


def main(fname_in, fname_out):
    vcf_reader = VCF(fname_in)
    vcf_writer = Writer(fname_out, vcf_reader)

    for variant in vcf_reader:
        variant.genotypes = [
            [*list(range(1, len(variant.ALT) + 1)), False]
        ]  # genotype 0 is reference

        vcf_writer.write_record(variant)

    vcf_writer.close()
    vcf_reader.close()


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
