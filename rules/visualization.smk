rule generate_web_visualization:
    input:
        consensus_file = "{dataset}/references/ref_majority.fasta",
        coverage_file = "variants/coverage.tsv",
        vcf_file = "{dataset}/variants/SNVs/snvs.vcf",
        gff_directory = config.input['gff_directory'] if config.input['gff_directory'] else [],
        primers_file = config.input['primers_file'] if config.input['primers_file'] else [],
        phylogeny_data = config.input['phylogeny_data'] if config.input['phylogeny_data'] else [],
        metainfo_file = config.input['metainfo_file'] if config.input['metainfo_file'] else [],
        global_ref = "variants/cohort_consensus.fasta", # see input.REF in rule snv 
        reference_file = "references/NC_045512.2.fasta",
        bam_file = "{dataset}/alignments/REF_aln.bam"
    output:
        alignment_file = "{dataset}/visualization/alignment.fasta",
        nwk_file = "{dataset}/visualization/tree.nwk",
        html_file = "{dataset}/visualization/index.html",
        reference_uri_file = "{dataset}/visualization/reference_uri_file",
        bam_uri_file = "{dataset}/visualization/bam_uri_file"
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

        augur align --sequences {input.phylogeny_data} {input.consensus_file} --output {output.alignment_file}
        augur tree --alignment {output.alignment_file} --output {output.nwk_file}

        create_datauri {input.reference_file} > {output.reference_uri_file}
        create_datauri {input.bam_file} > {output.bam_uri_file}

        python "{workflow.basedir}/scripts/assemble_visualization_webpage.py" \
            --consensus	"{input.consensus_file}" \
            --coverage	"{input.coverage_file}" \
            --vcf	"{input.vcf_file}" \
            --gff	"{input.gff_directory}" \
            --primers	"{input.primers_file}" \
            --nwk       "{output.nwk_file}" \
            --metainfo	"{input.metainfo_file}" \
            --template	"{workflow.basedir}/scripts/visualization.html" \
            --output	"{output.html_file}" \
            --wildcards	"{wildcards.dataset}" \
            --reference	"{input.global_ref}" \
            --reference_uri_file "{output.reference_uri_file}" \
            --bam_uri_file "{output.bam_uri_file}" \
            > {log.outfile} 2> {log.errfile}
        """
