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


rule generate_web_visualization:
    input:
        vcf_file = "{dataset}/variants/SNVs/snvs.vcf",
        gff_file = config.input['gff_file']
    output:
        html_file = "{dataset}/visualization/index.html"
    run:
        # load biodata in json format
        vcf_json = convert_vcf(input.vcf_file)
        gff_json = convert_gff(input.gff_file)

        embed_code = f"""
            var vcfData = {vcf_json}
            var gffData = {gff_json}
        """

        # assemble webpage
        raw_html_file = Path(workflow.basedir) / 'scripts' / 'visualization.html'
        with open(raw_html_file) as fd:
            raw_html = fd.read()

        # mod_html = raw_html.format(EXTERNAL_SNAKEMAKE_CODE_MARKER=embed_code)
        mod_html = raw_html.replace('{EXTERNAL_SNAKEMAKE_CODE_MARKER}', embed_code)

        with open(output.html_file, 'w') as fd:
            fd.write(mod_html)

    # seems to crash due to structure of config object: (TODO: investigate this)
    # script:
    #     '../scripts/assemble_visualization_webpage.py'
