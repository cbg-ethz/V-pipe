# NOTE  shorah and lofreq each have their independent option "consensus"
# see input.REF in rule snv
rule generate_web_visualization:
    input:
        consensus_file="{dataset}/references/ref_majority.fasta",
        coverage_file="{dataset}/alignments/coverage.tsv.gz",
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
        global_ref=(
            os.path.join(
                config.output["datadir"],
                config.output["cohortdir"],
                "cohort_consensus.fasta",
            )
            if config["lofreq" if config.general["snv_caller"] == "lofreq" else "snv"][
                "consensus"
            ]
            else reference_file
        ),
    output:
        html_file="{dataset}/visualization/index.html",
    params:
        assemble_visualization_webpage=srcdir(
            "../scripts/assemble_visualization_webpage.py"
        ),
        visualization_template=srcdir("../scripts/visualization.html"),
    log:
        outfile="{dataset}/visualization/stdout.log",
        errfile="{dataset}/visualization/stderr.log",
    conda:
        config.web_visualization["conda"]
    benchmark:
        "{dataset}/visualization/html_generation.benchmark"
    # group: 'snv' # HACK it's too fast and it confuses snakemake's timestamping
    resources:
        disk_mb=2000,
        mem_mb=config.web_visualization["mem"],
        time_min=config.web_visualization["time"],
    threads: 1
    shell:
        """
        # Why a shell directive?
        # 1) script directive crashes with `VpipeConfig`
        # 2) run directive does not allow conda envs

        python "{params.assemble_visualization_webpage}" \
            --consensus    "{input.consensus_file}" \
            --coverage    "{input.coverage_file}" \
            --vcf    "{input.vcf_file}" \
            --gff    "{input.gff_directory}" \
            --primers    "{input.primers_file}" \
            --metainfo    "{input.metainfo_file}" \
            --template    "{params.visualization_template}" \
            --output    "{output.html_file}" \
            --wildcards    "{wildcards.dataset}" \
            --reference    "{input.global_ref}" \
            > {log.outfile} 2> {log.errfile}
        """
