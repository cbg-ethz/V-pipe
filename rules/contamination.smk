import os

__author__ = "Susana Posada-Cespedes"
__author__ = "David Seifert"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"


# 1. align against 5VM as a QA check
rule bwa_QA:
    input:
        patient_ref="{dataset}/references/ref_{kind}.fasta",
        virusmix_ref=config.bwa_QA["ref_panel"],
        FASTQ=input_align,
    output:
        SAM=temp("{dataset}/QA_alignments/bwa_QA_{kind}.sam"),
        MSA="{dataset}/QA_alignments/bwa_refs_msa_{kind}.fasta",
    params:
        scratch="1250",
        mem=config.bwa_QA["mem"],
        time=config.bwa_QA["time"],
        BWA=config.applications["bwa"],
        MAFFT=config.applications["mafft"],
    log:
        outfile="{dataset}/QA_alignments/bwa_{kind}.out.log",
        errfile="{dataset}/QA_alignments/bwa_{kind}.err.log",
    conda:
        config.bwa_QA["conda"]
    benchmark:
        "{dataset}/QA_alignments/bwa_{kind}.benchmark"
    threads: config.bwa_QA["threads"]
    shell:
        """
        # 1. cleanup old run
        rm -f {output.SAM} {output.MSA}

        # 2. concatenate references
        mkdir -p {wildcards.dataset}/QA_alignments
        cat {input.patient_ref} {input.virusmix_ref} > {wildcards.dataset}/QA_alignments/bwa_refs_{wildcards.kind}.fasta

        # 3. indexing
        {params.BWA} index {wildcards.dataset}/QA_alignments/bwa_refs_{wildcards.kind}.fasta 2> >(tee {log.errfile} >&2)

        # 4. align
        {params.BWA} mem -t {threads} {wildcards.dataset}/QA_alignments/bwa_refs_{wildcards.kind}.fasta {input.FASTQ} > {output.SAM} 2> >(tee -a {log.errfile} >&2)

        # 5. MSA
        {params.MAFFT} --nuc --preservecase --maxiterate 1000 --localpair --thread {threads} {wildcards.dataset}/QA_alignments/bwa_refs_{wildcards.kind}.fasta > {output.MSA} 2> >(tee -a {log.errfile} >&2)

        # 6. cleanup BWA indices
        rm -f {wildcards.dataset}/QA_alignments/bwa_refs_{wildcards.kind}.fasta.*
        """


# 2. Call coverage statistics
rule coverage_QA:
    input:
        BAM="{dataset}/QA_alignments/bwa_QA_{kind}.bam",
        MSA="{dataset}/QA_alignments/bwa_refs_msa_{kind}.fasta",
    output:
        "{dataset}/QA_alignments/coverage_{kind}.tsv",
    params:
        scratch="1250",
        mem=config.coverage_QA["mem"],
        time=config.coverage_QA["time"],
        TARGET=config.coverage_QA["target"],
        COV_STATS=config.applications["coverage_stats"],
    log:
        outfile="{dataset}/QA_alignments/coverage_QA_{kind}.out.log",
        errfile="{dataset}/QA_alignments/coverage_QA_{kind}.err.log",
    conda:
        config.coverage_QA["conda"]
    benchmark:
        "{dataset}/QA_alignments/coverage_QA_{kind}.benchmark"
    threads: 1
    shell:
        """
        CONSENSUS_NAME={wildcards.dataset}
        CONSENSUS_NAME="${{CONSENSUS_NAME#*/}}"
        CONSENSUS_NAME="${{CONSENSUS_NAME//\//-}}"

        # 1. clean previous run
        rm -f {output}
        mkdir -p {wildcards.dataset}/QA_alignments

        # 2. collect coverage stats
        # we only collect statistics in the loop regions
        # of HIV-1 in order
        {params.COV_STATS} -t {params.TARGET} -i {input.BAM} -o {output} -m {input.MSA} --select "${{CONSENSUS_NAME}}" > {log.outfile} 2> >(tee {log.errfile} >&2)
        """
