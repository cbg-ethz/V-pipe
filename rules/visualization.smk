rule generate_web_visualization:
    input:
        consensus_file = "{dataset}/references/ref_majority.fasta",
        coverage_file = "variants/coverage.tsv",
        vcf_file = "{dataset}/variants/SNVs/snvs.vcf",
        gff_directory = config.input['gff_directory'],
        primers_file = config.input['primers_file'],
        global_ref = reference_file
    output:
        html_file = "{dataset}/visualization/index.html"
    params:
        scratch = "2000",
        mem = config.web_visualization["mem"],
        time = config.web_visualization["time"],
    log:
        outfile = "{dataset}/visualization/stdout.log",
        errfile = "{dataset}/visualization/stderr.log"
    conda:
        config.web_visualization["conda"]
    benchmark:
        "{dataset}/visualization/html_generation.benchmark"
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
            "{input.primers_file}" \
            "{workflow.basedir}/scripts/visualization.html" \
            "{output.html_file}" \
            "{wildcards.dataset}" \
            "{input.global_ref}" \
            > {log.outfile} 2> {log.errfile}
        """
