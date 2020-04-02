# TODO: paths probably need adjustment
rule generate_web_visualization:
    input:
        vcf_file = "{dataset}/variants/SNVs/snvs.vcf",  # NOTE: does not exist yet, will be merged in later
        gff_file = config.input['gff_file']
    output:
        html_file = "{dataset}/variants/visualization/index.html"  # NOTE: is not yet picked up by `rule all`
    script:
        'scripts/assemble_visualization_webpage.py'
