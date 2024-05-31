__author__ = "Kim"
__author__ = "Lara Fuhrmann"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"


rule consensus_bcftools:
    input:
        fname_bam=alignment_wildcard,
        fname_cov="{dataset}/alignments/coverage.tsv.gz",
        fname_ref=reference_file,
    output:
        fname_bcf="{dataset}/references/consensus.bcftools.bcf.gz",
        fname_temp_bcf=temp("{dataset}/references/temp.consensus.bcftools.bcf.gz"),
        fname_fasta="{dataset}/references/consensus.bcftools.fasta",
        fname_chain="{dataset}/references/consensus.bcftools.chain",
        fname_fasta_ambig="{dataset}/references/consensus_ambig.bcftools.fasta",
        fname_chain_ambig="{dataset}/references/consensus_ambig.bcftools.chain",
        fname_mask_lowcoverage=temp(
            "{dataset}/references/coverage_mask_lowcoverage.bed"
        ),
    params:
        tsvbased=config.general["tsvbased"],
        max_coverage=config.consensus_bcftools["max_coverage"],
        mask_coverage_threshold=config.consensus_bcftools["mask_coverage_threshold"],
        ambiguous_base_coverage_threshold=config.consensus_bcftools[
            "ambiguous_base_coverage_threshold"
        ],
        enhance_bcf=cachepath(
            "../scripts/enhance_bcf.py", executable=True, localsource=True
        ),
        gunzip=config.applications["gunzip"],
        bcftools=config.applications["bcftools"],
    log:
        outfile="{dataset}/references/consensus.bcftools.out.log",
        errfile="{dataset}/references/consensus.bcftools.err.log",
    conda:
        config.consensus_bcftools["conda"]
    benchmark:
        "{dataset}/alignments/consensus.bcftools.benchmark"
    resources:
        disk_mb=1250,
        mem_mb=config.consensus_bcftools["mem"],
        runtime=config.consensus_bcftools["time"],
    threads: config.consensus_bcftools["threads"]
    shell:
        """
        {params.bcftools} mpileup \
            --threads {threads} \
            -Ou \
            -f {input.fname_ref} \
            --max-depth {params.max_coverage} \
            --max-idepth {params.max_coverage} \
            --annotate FORMAT/AD,FORMAT/DP,INFO/AD \
            {input.fname_bam} \
        | {params.bcftools} call \
            --threads {threads} \
            -Ou \
            -mv \
            --keep-alts \
        | {params.bcftools} norm \
            --threads {threads} \
            -Ou \
            -f {input.fname_ref} \
        | {params.bcftools} filter \
            --threads {threads} \
            -e 'TYPE="INDEL" & INFO/AD[1]<INFO/AD[0]' \
            -Ob \
            --output {output.fname_temp_bcf} \
            2> >(tee {log.errfile} >&2)
        #bcftools csq -f {input.fname_ref} -g wheretogetthis.gff3.gz in.vcf -Ob -o out.bcf

        {params.gunzip} -c {input.fname_cov} | tail -n +2 \
        | awk -v base={params.tsvbased} \
            '$3 < {params.mask_coverage_threshold} {{printf "%s\\t%d\\t%d\\n", $1, $2 - base, $2 - base + 1}}' \
        > {output.fname_mask_lowcoverage} \
            2> >(tee -a {log.errfile} >&2)

        # preparations
        {params.enhance_bcf} \
           {output.fname_temp_bcf} \
           {output.fname_bcf} \
           {params.ambiguous_base_coverage_threshold} \
           2> >(tee -a {log.errfile} >&2)

        {params.bcftools} index {output.fname_bcf} 2> >(tee -a {log.errfile} >&2)

        common_consensus_params="--fasta-ref {input.fname_ref} --mark-del - --mask {output.fname_mask_lowcoverage} --mask-with n"

        # majority bases
        {params.bcftools} consensus \
            --output {output.fname_fasta} \
            --chain {output.fname_chain} \
            $common_consensus_params \
            -H A \
            -i "INFO/AD[0]<INFO/AD[*]" \
            {output.fname_bcf} \
            2> >(tee -a {log.errfile} >&2)

        # ambiguous bases
        {params.bcftools} consensus \
            --output {output.fname_fasta_ambig} \
            --chain {output.fname_chain_ambig} \
            $common_consensus_params \
            -H I --iupac-codes \
            {output.fname_bcf} \
            2> >(tee -a {log.errfile} >&2)
        """


rule consensus_sequences:
    input:
        BAM=alignment_wildcard,
        REF=reference_file,
    output:
        REF_amb="{dataset}/references/ref_ambig.fasta",
        REF_amb_dels="{dataset}/references/ref_ambig_dels.fasta",
        REF_majority="{dataset}/references/ref_majority.fasta",
        REF_majority_dels="{dataset}/references/ref_majority_dels.fasta",
    params:
        MIN_COVERAGE=config.consensus_sequences["min_coverage"],
        N_COVERAGE=config.consensus_sequences["n_coverage"],
        QUAL_THRD=config.consensus_sequences["qual_thrd"],
        MIN_FREQ=config.consensus_sequences["min_freq"],
        OUTDIR="{dataset}/references",
        EXTRACT_CONSENSUS=config.applications["extract_consensus"],
    log:
        outfile="{dataset}/references/consensus_sequences.out.log",
        errfile="{dataset}/references/consensus_sequences.err.log",
    conda:
        config.consensus_sequences["conda"]
    benchmark:
        "{dataset}/alignments/consensus.benchmark"
    resources:
        disk_mb=1250,
        mem_mb=config.consensus_sequences["mem"],
        runtime=config.consensus_sequences["time"],
    threads: 1
    shell:
        """
        CONSENSUS_NAME={wildcards.dataset}
        CONSENSUS_NAME="${{CONSENSUS_NAME#*/}}"
        CONSENSUS_NAME="${{CONSENSUS_NAME//\//-}}"

        {params.EXTRACT_CONSENSUS} -i {input.BAM} -f {input.REF} -c {params.MIN_COVERAGE} -n {params.N_COVERAGE} -q {params.QUAL_THRD} -a {params.MIN_FREQ} -N "${{CONSENSUS_NAME}}" -o {params.OUTDIR} > {log.outfile} 2> >(tee -a {log.errfile} >&2)
        """


# QA checks performed on consensus_sequences
# - do pairwise alignement
rule consseq_QA:
    input:
        REF=reference_file,
        REF_majority_dels="{dataset}/references/ref_majority_dels.fasta",
    output:
        REF_matcher="{dataset}/references/ref_majority_dels.matcher",
    params:
        MATCHER=config.applications["matcher"],
    log:
        outfile="{dataset}/references/qa_consseq.out.log",
        errfile="{dataset}/references/qa_consseq.err.log",
    conda:
        config.consseq_QA["conda"]
    benchmark:
        "{dataset}/alignments/qa_consseq.benchmark"
    resources:
        disk_mb=1250,
        mem_mb=config.consseq_QA["mem"],
        runtime=config.consseq_QA["time"],
    threads: 1
    shell:
        """
        if tail -n +2 {input.REF_majority_dels} | grep -qE '[^n]'; then
            {params.MATCHER} -asequence {input.REF} -bsequence {input.REF_majority_dels} -outfile {output.REF_matcher} 2> >(tee {log.errfile} >&2)
        else
            touch {output.REF_matcher}
            echo "pure 'nnnn...' consensus, no possible alignement" | tee {log.outfile}
        fi
        """


rule frameshift_deletions_checks:
    input:
        REF_NAME=reference_file,
        BAM=alignment_wildcard,
        CONSENSUS="{dataset}/references/consensus.bcftools.fasta",
        CHAIN="{dataset}/references/consensus.bcftools.chain",
        GENES_GFF=(
            cachepath(config.input["genes_gff"]) if config.input["genes_gff"] else []
        ),
    output:
        FRAMESHIFT_DEL_CHECK_TSV="{dataset}/references/frameshift_deletions_check.tsv",
    params:
        FRAMESHIFT_DEL_CHECKS=config.applications["frameshift_deletions_checks"],
    log:
        outfile="{dataset}/references/frameshift_deletions_check.out.log",
        errfile="{dataset}/references/frameshift_deletions_check.err.log",
    conda:
        config.frameshift_deletions_checks["conda"]
    benchmark:
        "{dataset}/alignments/frameshift_deletions_check.benchmark"
    resources:
        disk_mb=1250,
        mem_mb=config.frameshift_deletions_checks["mem"],
        runtime=config.frameshift_deletions_checks["time"],
    threads: 1
    shell:
        """
        {params.FRAMESHIFT_DEL_CHECKS} -i {input.BAM} -c {input.CONSENSUS} --chain {input.CHAIN} -f {input.REF_NAME} -g {input.GENES_GFF} --english=true -o {output.FRAMESHIFT_DEL_CHECK_TSV} 2> >(tee {log.errfile} >&2)
        """


if config.general["aligner"] == "ngshmmalign":

    ruleorder: hmm_align > consensus_sequences


#
#
# elif config.general["aligner"] == "bwa":
#
#    ruleorder: consensus_sequences > hmm_align
#
#
# elif config.general["aligner"] == "bowtie":
#
#    ruleorder: consensus_sequences > hmm_align
