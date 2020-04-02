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


def main(vcf_file, gff_file, html_file):
    # load biodata in json format
    vcf_json = convert_vcf(vcf_file)
    gff_json = convert_gff(gff_file)

    embed_code = f"""
        var vcfData = {vcf_json}
        var gffData = {gff_json}
    """

    # assemble webpage
    raw_html_file = Path(workflow.basedir) / 'scripts' / 'visualization.html'
    with open(raw_html_file) as fd:
        raw_html = fd.read()

    raw_html.format(EXTERNAL_SNAKEMAKE_CODE_MARKER=embed_code)

    with open(html_file, 'w') as fd:
        fd.write(raw_html)


if __name__ == '__main__':
    main(
        snakemake.input.vcf_file, snakemake.input.gff_file,
        snakemake.output.html_file
    )
