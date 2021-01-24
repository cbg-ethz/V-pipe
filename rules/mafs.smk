
__author__ = "Susana Posada-Cespedes"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"

# TODO move and merge with consensus_sequences in align.smk
# 0. Gather base counts from alignments
rule basecounts:
    input:
        BAM = "{dataset}/alignments/REF_aln.bam",
    output:
        BASECNT = "{dataset}/alignments/basecnt.tsv.gz",
        COVERAGE = "{dataset}/alignments/coverage.tsv.gz",
        STATS = "{dataset}/alignments/REF_aln_stats.yaml",
    params:
        NAME = ID,
        ALN2BASECNT = config.applications['aln2basecnt'],
    log:
        outfile = "{dataset}/alignments/basecounts.out.log",
        errfile = "{dataset}/alignments/basecounts.out.log",
    conda:
        config.basecounts['conda']
    benchmark:
        "{dataset}/alignments/basecounts.benchmark"
    resources:
        disk_mb = 1250,
        mem_mb = config.basecounts['mem'],
        time_min = config.basecounts['time'],
    threads:
        1
    shell:
        """
        {params.ALN2BASECNT} --basecnt "{output.BASECNT}" --coverage "{output.COVERAGE}" --name "{params.NAME}" --stats "{output.STATS}" "{input.BAM}"
        """

# 1. Gather coverages into central big file
localrules:
    coverage_list
rule coverage_list:
    input:
        expand("{dataset}/alignments/coverage.tsv.gz", dataset=datasets),
    output:
        temp("variants/coverage.tmp.list"),
    run:
        with open(output[0], 'w') as out:
            out.write('\n'.join(input))

#HACK: troubles passing command lines to bash using `-c` that exceed 128kB (MAX_ARGS in limits.h) even when `getconf ARG_MAX` gives higher limits. We use a @list file instead.
rule coverage:
    input:
        SAMPLECOVS = expand("{dataset}/alignments/coverage.tsv.gz", dataset=datasets),
        COVLIST = ancient("variants/coverage.tmp.list"),
    output:
        COVERAGE = "variants/coverage.tsv",
        COVSTATS = "variants/coverage_stats.tsv",
    params:
        GATHER_COVERAGE = config.applications['gather_coverage'],
    log:
        outfile = "variants/coverage.out.log",
        errfile = "variants/coverage.out.log",
    conda:
        config.coverage['conda']
    benchmark:
        "variants/minority_variants.benchmark"
    resources:
        disk_mb = 1250,
        mem_mb = config.coverage['mem'],
        time_min = config.coverage['time'],
    threads:
        config.coverage['threads']
    shell:
        """
        {params.GATHER_COVERAGE} --output {output.COVERAGE} --stats {output.COVSTATS} --threads {threads} @{input.COVLIST} > >(tee {log.outfile}) 2>&1
        """

# TODO rewrite/adapt
# 2. Ouptut minor allele frequencies
rule minor_variants:
    input:
        REF = reference_file,
        BASECNT = expand("{dataset}/alignments/basecnt.tsv.gz", dataset=datasets),
    output:
        VARIANTS = "variants/minority_variants.tsv",
        CONSENSUS = "variants/cohort_consensus.fasta",
    params:
        OUTDIR = "variants",
        NAMES = IDs,
        MIN_COVERAGE = config.minor_variants['min_coverage'],
        FREQUENCIES = '--freqs' if config.minor_variants['frequencies'] else '',
        MINORITY_CALLER = config.applications['minority_freq'],
    log:
        outfile = "variants/minority_variants.out.log",
        errfile = "variants/minority_variants.out.log",
    conda:
        config.minor_variants['conda']
    benchmark:
        "variants/minority_variants.benchmark"
    resources:
        disk_mb = 1250,
        mem_mb = config.minor_variants['mem'],
        time_min = config.minor_variants['time'],
    threads:
        config.minor_variants['threads']
    shell:
        """
        exit 2
        """
        #{params.MINORITY_CALLER} -r {input.REF} -c {params.MIN_COVERAGE} -N {params.NAMES} -t {threads} -o {params.OUTDIR} -d {params.FREQUENCIES} {input.BAM} > >(tee {log.outfile}) 2>&1
