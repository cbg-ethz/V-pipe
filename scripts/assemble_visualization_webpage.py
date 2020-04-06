import sys
import json
from pathlib import Path

import vcf
from BCBio import GFF


def convert_vcf(fname):
    """Convert VCF to JSON."""
    output = []

    with open(fname) as fd:
        vcf_reader = vcf.Reader(fd)

        for record in vcf_reader:
            output.append({
                'position': record.POS,
                'variant': [v.sequence for v in record.ALT],
                'pvalue': record.INFO['Pval']
            })

    return json.dumps(output)


def convert_gff(fname):
    """Convert GFF to JSON."""
    output = {
        'genomeLength': None,
        'features': []
    }

    with open(fname) as fd:
        for record in GFF.parse(fd):
            assert output['genomeLength'] is None
            output['genomeLength'] = len(record.seq)

            for feature in record.features:
                output['features'].append({
                    'id': record.id,
                    'type': feature.type,
                    'name': feature.id,
                    'start': int(feature.location.start),
                    'end': int(feature.location.end)
                })

    return json.dumps(output)


def main(vcf_file, gff_file, html_file_in, html_file_out):
    # load biodata in json format
    vcf_json = convert_vcf(vcf_file)
    gff_json = convert_gff(gff_file)

    embed_code = f"""
        var vcfData = {vcf_json}
        var gffData = {gff_json}
    """

    # assemble webpage
    with open(html_file_in) as fd:
        raw_html = fd.read()

    # TODO: make this more robust
    mod_html = raw_html.replace('{EXTERNAL_SNAKEMAKE_CODE_MARKER}', embed_code)

    with open(html_file_out, 'w') as fd:
        fd.write(mod_html)


if __name__ == '__main__':
    main(
        sys.argv[1], sys.argv[2],
        sys.argv[3], sys.argv[4]
    )
