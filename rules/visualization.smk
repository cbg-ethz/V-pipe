rule generate_web_visualization:
    input:
        vcf_file = "{dataset}/variants/SNVs/snvs.vcf",
        gff_file = config.input['gff_file']
    output:
        html_file = "{dataset}/visualization/index.html"
    conda:
        '../envs/visualization.yaml'
    shell:
        """
        # Why a shell directive?
        # 1) script directive crashes with `VpipeConfig`
        # 2) run directive does not allow conda envs

        python "{workflow.basedir}/scripts/assemble_visualization_webpage.py" \
            "{input.vcf_file}" \
            "{input.gff_file}" \
            "{workflow.basedir}/scripts/visualization.html" \
            "{output.html_file}"
        """
