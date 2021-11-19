# NOTE  shorah and lofreq each have their independent option "consensus"
# see input.REF in rule snv
rule generate_web_visualization:
    input:
        consensus_file="{dataset}/references/ref_majority.fasta",
        reference_file = config["input"]["reference"],
        coverage_file="{dataset}/alignments/coverage.tsv.gz",
        vcf_file="{dataset}/variants/SNVs/snvs.vcf",
        bam_file="{dataset}/alignments/REF_aln.bam",
        gff_directory=(
            config.input["gff_directory"] if config.input["gff_directory"] else []
        ),
        primers_file=(
            config.input["primers_file"] if config.input["primers_file"] else []
        ),
        metainfo_file=(
            config.input["metainfo_file"] if config.input["metainfo_file"] else []
        ),
        phylogeny_data=(
            config.input["phylogeny_data"] if config.input["phylogeny_data"] else []
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
        snv_html_file="{dataset}/visualization/snv_calling.html",
        alignment_html_file="{dataset}/visualization/alignment.html",
        reference_uri_file="{dataset}/visualization/reference_uri_file",
        bam_uri_file="{dataset}/visualization/bam_uri_file",
    params:
        assemble_visualization_webpage=cachepath(
            "../scripts/assemble_visualization_webpage.py",
            executable=True,
            localsource=True,
        ),
        snv_visualization_template=cachepath(
            "../scripts/snv_calling_visualization.html", localsource=True
        ),
        alignment_visualization_template=cachepath(
            "../scripts/alignment_visualization.html", localsource=True
        ),
        plot_phylogenetic_tree=config.input["phylogeny_data"],
        alignment_file=("{dataset}/visualization/alignment.fasta" if config.input["phylogeny_data"] else []),
        nwk_file=("{dataset}/visualization/tree.nwk" if config.input["phylogeny_data"] else []),
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

        if [ ! -z "{params.plot_phylogenetic_tree}" ]
        then
          # generate phylogenetic tree
          augur align --sequences {input.phylogeny_data} {input.consensus_file} --output {params.alignment_file}
          augur tree --alignment {params.alignment_file} --output {params.nwk_file}
        fi

        # generate read alignment
        create_datauri {input.reference_file} > {output.reference_uri_file}
        create_datauri {input.bam_file} > {output.bam_uri_file}

        # generate html visualization pages
        python "{params.assemble_visualization_webpage}" \
            --consensus    "{input.consensus_file}" \
            --coverage    "{input.coverage_file}" \
            --vcf    "{input.vcf_file}" \
            --gff    "{input.gff_directory}" \
            --primers    "{input.primers_file}" \
            --metainfo    "{input.metainfo_file}" \
            --snv_calling_template    "{params.snv_visualization_template}" \
            --alignment_template    "{params.alignment_visualization_template}" \
            --html_out_snv_calling    "{output.snv_html_file}" \
            --html_out_alignment    "{output.alignment_html_file}" \
            --wildcards    "{wildcards.dataset}" \
            --reference    "{input.global_ref}" \
            --reference_uri_file "{output.reference_uri_file}" \
            --bam_uri_file "{output.bam_uri_file}" \
            --nwk "{params.nwk_file}" \
            > {log.outfile} 2> {log.errfile}
       
        """
