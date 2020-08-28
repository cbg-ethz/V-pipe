__author__ = "Susana Posada-Cespedes"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"


# TODO move and merge with consensus_sequences in align.smk
# 0. Gather base counts from alignments
rule basecounts:
    input:
        BAM="{dataset}/alignments/REF_aln.bam",
    output:
        BASECNT="{dataset}/alignments/basecnt.tsv.gz",
        COVERAGE="{dataset}/alignments/coverage.tsv.gz",
        STATS="{dataset}/alignments/REF_aln_stats.yaml",
    params:
        scratch="1250",
        mem=config.basecounts["mem"],
        time=config.basecounts["time"],
        NAME=ID,
        ALN2BASECNT=config.applications["aln2basecnt"],
    log:
        outfile="{dataset}/alignments/basecounts.out.log",
        errfile="{dataset}/alignments/basecounts.out.log",
    conda:
        config.basecounts["conda"]
    benchmark:
        "{dataset}/alignments/basecounts.benchmark"
    threads: 1
    shell:
        """
        {params.ALN2BASECNT} --basecnt "{output.BASECNT}" --coverage "{output.COVERAGE}" --name "{params.NAME}" --stats "{output.STATS}" "{input.BAM}"
        """


# 1. Gather coverages into central big file
rule coverage:
    input:
        SAMPLECOVS=expand("{dataset}/alignments/coverage.tsv.gz", dataset=datasets),
    output:
        COVERAGE="variants/coverage.tsv",
        COVSTATS="variants/coverage_stats.tsv",
    params:
        scratch="1250",
        mem=config.coverage["mem"],
        time=config.coverage["time"],
        GATHER_COVERAGE=config.applications["gather_coverage"],
    log:
        outfile="variants/coverage.out.log",
        errfile="variants/coverage.out.log",
    conda:
        config.coverage["conda"]
    benchmark:
        "variants/minority_variants.benchmark"
    threads: config.coverage["threads"]
    shell:
        """
        {params.GATHER_COVERAGE} --output {output.COVERAGE} --stats {output.COVSTATS} --threads {threads} {input.SAMPLECOVS} > >(tee {log.outfile}) 2>&1
        """


# TODO rewrite/adapt
# 2. Ouptut minor allele frequencies
rule minor_variants:
    input:
        REF=reference_file,
        BASECNT=expand("{dataset}/alignments/basecnt.tsv.gz", dataset=datasets),
    output:
        VARIANTS="variants/minority_variants.tsv",
        CONSENSUS="variants/cohort_consensus.fasta",
    params:
        scratch="1250",
        mem=config.minor_variants["mem"],
        time=config.minor_variants["time"],
        OUTDIR="variants",
        NAMES=IDs,
        MIN_COVERAGE=config.minor_variants["min_coverage"],
        FREQUENCIES="--freqs" if config.minor_variants["frequencies"] else "",
        MINORITY_CALLER=config.applications["minority_freq"],
    log:
        outfile="variants/minority_variants.out.log",
        errfile="variants/minority_variants.out.log",
    conda:
        config.minor_variants["conda"]
    benchmark:
        "variants/minority_variants.benchmark"
    threads: config.minor_variants["threads"]
    shell:
        """
        exit 2
        """
        #{params.MINORITY_CALLER} -r {input.REF} -c {params.MIN_COVERAGE} -N {params.NAMES} -t {threads} -o {params.OUTDIR} -d {params.FREQUENCIES} {input.BAM} > >(tee {log.outfile}) 2>&1
