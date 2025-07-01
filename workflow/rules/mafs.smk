import os

__author__ = "Susana Posada-Cespedes"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"


# TODO move and merge with consensus_sequences in align.smk
# 0. Gather base counts from alignments
rule basecounts:
    input:
        BAM=alignment_wildcard,
    output:
        BASECNT="{dataset}/alignments/basecnt.tsv.gz",
        COVERAGE="{dataset}/alignments/coverage.tsv.gz",
        STATS="{dataset}/alignments/REF_aln_stats.yaml",
    params:
        NAME=ID,
        ALN2BASECNT=config.applications["aln2basecnt"],
        ARRAYBASED=config.general["tsvbased"],
    log:
        outfile="{dataset}/alignments/basecounts.out.log",
        errfile="{dataset}/alignments/basecounts.out.log",
    conda:
        config.basecounts["conda"]
    benchmark:
        "{dataset}/alignments/basecounts.benchmark"
    resources:
        disk_mb=1250,
        mem_mb=config.basecounts["mem"],
        runtime=config.basecounts["time"],
    threads: 1
    shell:
        """
        {params.ALN2BASECNT} --first "{params.ARRAYBASED}" --basecnt "{output.BASECNT}" --coverage "{output.COVERAGE}" --name "{params.NAME}" --stats "{output.STATS}" "{input.BAM}" > {log.outfile} 2> >(tee {log.errfile} >&2)
        """


rule chromsize:
    input:
        reference_file,
    output:
        chrom_size=cohortdir("chrom.size"),
    params:
        CHROMSIZE=config.applications["chromsize"],
    log:
        outfile=cohortdir("chromsize.out.log"),
        errfile=cohortdir("chromsize.err.log"),
    conda:
        config.chromsize["conda"]
    benchmark:
        cohortdir("chromsize.benchmark")
    resources:
        disk_mb=1250,
        mem_mb=config.chromsize["mem"],
        runtime=config.chromsize["time"],
    threads: 1
    shell:
        r"""
        {params.CHROMSIZE} --accession-only --fasta "{input}" --output "{output.chrom_size}" \
            > {log.outfile} 2> >(tee -a "{log.errfile}" >&2)
        """


rule basecounts_QC:
    input:
        COVERAGE="{dataset}/alignments/coverage.tsv.gz",
        CHROM_SIZE=(
            cohortdir("chrom.size")
            if config["basecounts_qc"]["depth_qc_type"] == "fraction"
            else []
        ),
    output:
        COV_DEPTH_QC="{dataset}/alignments/coverage_depth_qc.yaml",
    params:
        COV_DEPTH_QC=config.applications["coverage_depth_qc"],
        DEPTHS=config["basecounts_qc"]["depth_qc_list"],
        CHROM_SIZE=(
            f"--fract={cohortdir('chrom.size')}"
            if config["basecounts_qc"]["depth_qc_type"] == "fraction"
            else ""
        ),
    log:
        outfile="{dataset}/alignments/basecounts_qc.out.log",
        errfile="{dataset}/alignments/basecounts_qc.out.log",
    conda:
        config.basecounts_qc["conda"]
    benchmark:
        "{dataset}/alignments/coverage_depth_qc.benchmark"
    resources:
        disk_mb=1250,
        mem_mb=config.basecounts_qc["mem"],
        runtime=config.basecounts_qc["time"],
    threads: 1
    shell:
        """
        {params.COV_DEPTH_QC} {params.CHROM_SIZE} --depth {params.DEPTHS} --output {output.COV_DEPTH_QC} -- {input.COVERAGE}    \
            > "{log.outfile}" 2> >(tee "{log.errfile}" >&2)
        """


# 1. Gather coverages into central big file
localrules:
    coverage_list,


rule coverage_list:
    input:
        expand("{dataset}/alignments/coverage.tsv.gz", dataset=datasets),
    output:
        temp(cohortdir("coverage.tmp.list")),
    run:
        with open(output[0], "w") as out:
            out.write("\n".join(input))


# HACK: troubles passing command lines to bash using `-c` that exceed 128kB (MAX_ARGS in limits.h) even when `getconf ARG_MAX` gives higher limits. We use a @list file instead.
rule coverage:
    input:
        SAMPLECOVS=expand("{dataset}/alignments/coverage.tsv.gz", dataset=datasets),
        COVLIST=ancient(cohortdir("coverage.tmp.list")),
    output:
        COVERAGE=cohortdir("coverage.tsv"),
        COVSTATS=cohortdir("coverage_stats.tsv"),
    params:
        GATHER_COVERAGE=config.applications["gather_coverage"],
    log:
        outfile=cohortdir("coverage.out.log"),
        errfile=cohortdir("coverage.out.log"),
    conda:
        config.coverage["conda"]
    benchmark:
        cohortdir("minority_variants.benchmark")
    resources:
        disk_mb=1250,
        mem_mb=config.coverage["mem"],
        runtime=config.coverage["time"],
    threads: config.coverage["threads"]
    shell:
        """
        {params.GATHER_COVERAGE} --output {output.COVERAGE} --stats {output.COVSTATS} --threads {threads} @{input.COVLIST} > >(tee {log.outfile}) 2> >(tee {log.errfile} >&2)
        """


# TODO rewrite/adapt
# 2. Ouptut minor allele frequencies
rule minor_variants:
    input:
        REF=reference_file,
        BASECNT=expand("{dataset}/alignments/basecnt.tsv.gz", dataset=datasets),
        BAM=expand(alignment_wildcard, dataset=datasets),
    output:
        VARIANTS=cohortdir("minority_variants.tsv"),
        CONSENSUS=cohortdir("cohort_consensus.fasta"),
    params:
        OUTDIR=cohortdir(""),
        NAMES=IDs,
        MIN_COVERAGE=config.minor_variants["min_coverage"],
        FREQUENCIES="--freqs" if config.minor_variants["frequencies"] else "",
        MINORITY_CALLER=config.applications["minority_freq"],
    log:
        outfile=cohortdir("minority_variants.out.log"),
        errfile=cohortdir("minority_variants.out.log"),
    conda:
        config.minor_variants["conda"]
    benchmark:
        cohortdir("minority_variants.benchmark")
    resources:
        disk_mb=1250,
        mem_mb=config.minor_variants["mem"],
        runtime=config.minor_variants["time"],
    threads: config.minor_variants["threads"]
    shell:
        """
        {params.MINORITY_CALLER} -r {input.REF} -c {params.MIN_COVERAGE} -N {params.NAMES} -t {threads} -o {params.OUTDIR} {params.FREQUENCIES} {input.BAM} > >(tee {log.outfile}) 2> >(tee {log.errfile} >&2)
        """
