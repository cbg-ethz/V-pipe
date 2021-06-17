rule generate_web_visualization:
    input:
        consensus_file="{dataset}/references/ref_majority.fasta",
        coverage_file="variants/coverage.tsv",
        vcf_file="{dataset}/variants/SNVs/snvs.vcf",
        gff_directory=(
            config.input["gff_directory"] if config.input["gff_directory"] else []
        ),
        primers_file=(
            config.input["primers_file"] if config.input["primers_file"] else []
        ),
        metainfo_file=(
            config.input["metainfo_file"] if config.input["metainfo_file"] else []
        ),
        # see input.REF in rule snv
        global_ref=(
            "variants/cohort_consensus.fasta"
            if config.snv["consensus"]
            else reference_file
        ),
    output:
        html_file="{dataset}/visualization/index.html",
    params:
        scratch="2000",
        mem=config.web_visualization["mem"],
        time=config.web_visualization["time"],
    log:
        outfile="{dataset}/visualization/stdout.log",
        errfile="{dataset}/visualization/stderr.log",
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
            --consensus    "{input.consensus_file}" \
            --coverage    "{input.coverage_file}" \
            --vcf    "{input.vcf_file}" \
            --gff    "{input.gff_directory}" \
            --primers    "{input.primers_file}" \
            --metainfo    "{input.metainfo_file}" \
            --template    "{workflow.basedir}/scripts/visualization.html" \
            --output    "{output.html_file}" \
            --wildcards    "{wildcards.dataset}" \
            --reference    "{input.global_ref}" \
            > {log.outfile} 2> {log.errfile}
        """
