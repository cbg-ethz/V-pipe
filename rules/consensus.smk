rule consensus_bcftools:
    input:
        fname_bam="{dataset}/alignments/REF_aln.bam",
        fname_cov="{dataset}/alignments/coverage.tsv.gz",
        fname_ref=reference_file,
    output:
        fname_bcf="{dataset}/references/consensus.bcftools.bcf.gz",
        fname_temp_bcf=temp("{dataset}/references/temp.consensus.bcftools.bcf.gz"),
        fname_fasta="{dataset}/references/consensus.bcftools.fasta",
        fname_fasta_ambig="{dataset}/references/consensus_ambig.bcftools.fasta",
        fname_mask_lowcoverage=temp(
            "{dataset}/references/coverage_mask_lowcoverage.bed"
        ),
    params:
        max_coverage=config.consensus_bcftools["max_coverage"],
        mask_coverage_threshold=config.consensus_bcftools["mask_coverage_threshold"],
        ambiguous_base_coverage_threshold=config.consensus_bcftools[
            "ambiguous_base_coverage_threshold"
        ],
        script_dir=os.path.join(VPIPE_BASEDIR, "scripts"),
    log:
        outfile="{dataset}/references/consensus.bcftools.out.log",
        errfile="{dataset}/references/consensus.bcftools.err.log",
    conda:
        config.consensus_bcftools["conda"]
    resources:
        disk_mb=1250,
        mem_mb=config.consensus_bcftools["mem"],
        time_min=config.consensus_bcftools["time"],
    threads: config.consensus_bcftools["threads"]
    shell:
        """
        bcftools mpileup \
            --threads {threads} \
            -Ou \
            -f {input.fname_ref} \
            --max-depth {params.max_coverage} \
            --max-idepth {params.max_coverage} \
            --annotate FORMAT/AD,FORMAT/DP,INFO/AD \
            {input.fname_bam} \
        | bcftools call \
            --threads {threads} \
            -Ou \
            -mv \
            --keep-alts \
        | bcftools norm \
            --threads {threads} \
            -Ou \
            -f {input.fname_ref} \
        | bcftools filter \
            --threads {threads} \
            -e 'TYPE="INDEL" & INFO/AD[1]<INFO/AD[0]' \
            -Ob \
            --output {output.fname_temp_bcf}
        #bcftools csq -f {input.fname_ref} -g wheretogetthis.gff3.gz in.vcf -Ob -o out.bcf

        # TODO: homogene use of 0-base vs 1-base
        zcat {input.fname_cov} | tail -n +2 \
        | awk -v base=0 \
            '$3 < {params.mask_coverage_threshold} {{printf "%s\\t%d\\t%d\\n", $1, $2 - base, $2 - base + 1}}' \
        > {output.fname_mask_lowcoverage}

        preparations
        {params.script_dir}/enhance_bcf.py \
           {output.fname_temp_bcf} \
           {output.fname_bcf} \
           {params.ambiguous_base_coverage_threshold}

        bcftools index {output.fname_bcf}

        common_consensus_params="--fasta-ref {input.fname_ref} --mark-del - --mask {output.fname_mask_lowcoverage} --mask-with n"

        # majority bases
        bcftools consensus \
            --output {output.fname_fasta} \
            $common_consensus_params \
            -H A \
            -i "INFO/AD[0]<INFO/AD[*]" \
            {output.fname_bcf}

        # ambiguous bases
        bcftools consensus \
            --output {output.fname_fasta_ambig} \
            $common_consensus_params \
            -H I --iupac-codes \
            {output.fname_bcf}
        """
