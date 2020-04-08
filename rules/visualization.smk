rule generate_web_visualization:
    input:
        consensus_file = "{dataset}/references/ref_majority.fasta",
        coverage_file = "{dataset}/variants/coverage.tsv"
        vcf_file = "{dataset}/variants/SNVs/snvs.vcf",
        gff_directory = config.input['gff_directory']
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
            "{input.consensus_file}" \
            "{input.coverage_file}" \
            "{input.vcf_file}" \
            "{input.gff_directory}" \
            "{workflow.basedir}/scripts/visualization.html" \
            "{output.html_file}"
        """
