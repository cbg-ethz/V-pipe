rule consensus_bcftools:
    input:
        fname_bam = "{dataset}/alignments/REF_aln.bam",
        fname_ref = reference_file,
    output:
        fname_bcf = "{dataset}/references/consensus.bcftools.bcf.gz",
        fname_fasta = "{dataset}/references/consensus.bcftools.fasta",
        fname_fasta_ambig = "{dataset}/references/consensus_ambig.bcftools.fasta",
        fname_mask_lowcoverage = temp("{dataset}/references/coverage_mask_lowcoverage.bed"),
    params:
        scratch = '1250',
        mem = config.consensus_bcftools['mem'],
        time = config.consensus_bcftools['time'],

        max_coverage = config.consensus_bcftools['max_coverage'],
        mask_coverage_threshold = config.consensus_bcftools['mask_coverage_threshold'],
        ambiguous_base_coverage_threshold = config.consensus_bcftools['ambiguous_base_coverage_threshold'],

        script_dir = os.path.join(VPIPE_BASEDIR, 'scripts'),
    conda:
        config.consensus_bcftools['conda']
    threads: config.consensus_bcftools['threads']
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
            --output {output.fname_bcf}
        #bcftools csq -f {input.fname_ref} -g wheretogetthis.gff3.gz in.vcf -Ob -o out.bcf

        # TODO: retrieve this data from other rules
        common_samtools_params="-a -d 0 -J {input.fname_bam}"
        samtools depth \
            $common_samtools_params \
        | awk \
            '$3 < {params.mask_coverage_threshold} {{printf "%s\\t%d\\t%d\\n", $1, $2 - 1, $2}}' \
        > {output.fname_mask_lowcoverage}

        # preparations
        python3 {params.script_dir}/enhance_bcf.py \
            {output.fname_bcf} \
            temp.bcf.gz \
            {params.ambiguous_base_coverage_threshold}
        mv temp.bcf.gz {output.fname_bcf}

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
