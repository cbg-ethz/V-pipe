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
    log:
        outfile="{dataset}/alignments/basecounts.out.log",
        errfile="{dataset}/alignments/basecounts.out.log",
    benchmark:
        "{dataset}/alignments/basecounts.benchmark"
    conda:
        config.basecounts["conda"]
    threads: 1
    resources:
        disk_mb=1250,
        mem_mb=config.basecounts["mem"],
        runtime=config.basecounts["time"],
    params:
        NAME=ID,
        ALN2BASECNT=config.applications["aln2basecnt"],
        ARRAYBASED=config.general["tsvbased"],
    shell:
        """
        {params.ALN2BASECNT} --first "{params.ARRAYBASED}" --basecnt "{output.BASECNT}" --coverage "{output.COVERAGE}" --name "{params.NAME}" --stats "{output.STATS}" "{input.BAM}" > {log.outfile} 2> >(tee {log.errfile} >&2)
        """


rule chromsize:
    input:
        reference_file,
    output:
        chrom_size=cohortdir("chrom.size"),
    log:
        outfile=cohortdir("chromsize.out.log"),
        errfile=cohortdir("chromsize.err.log"),
    benchmark:
        cohortdir("chromsize.benchmark")
    conda:
        config.chromsize["conda"]
    threads: 1
    resources:
        disk_mb=1250,
        mem_mb=config.chromsize["mem"],
        runtime=config.chromsize["time"],
    params:
        CHROMSIZE=config.applications["chromsize"],
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
    log:
        outfile="{dataset}/alignments/basecounts_qc.out.log",
        errfile="{dataset}/alignments/basecounts_qc.out.log",
    benchmark:
        "{dataset}/alignments/coverage_depth_qc.benchmark"
    conda:
        config.basecounts_qc["conda"]
    threads: 1
    resources:
        disk_mb=1250,
        mem_mb=config.basecounts_qc["mem"],
        runtime=config.basecounts_qc["time"],
    params:
        COV_DEPTH_QC=config.applications["coverage_depth_qc"],
        DEPTHS=config["basecounts_qc"]["depth_qc_list"],
        CHROM_SIZE=(
            f"--fract={cohortdir('chrom.size')}"
            if config["basecounts_qc"]["depth_qc_type"] == "fraction"
            else ""
        ),
    shell:
        """
        {params.COV_DEPTH_QC} {params.CHROM_SIZE} --depth {params.DEPTHS} --output {output.COV_DEPTH_QC} -- {input.COVERAGE}    \
            > "{log.outfile}" 2> >(tee "{log.errfile}" >&2)
        """


rule classif_by_coverage:
    input:
        COVERAGE="{dataset}/alignments/coverage.tsv.gz",
        CLASSIF_BED=config["classif_by_coverage"]["bed_file"],
    output:
        CLASSIF_CSV="{dataset}/alignments/classif_by_coverage.csv",
    log:
        outfile="{dataset}/alignments/classif_by_coverage.out.log",
        errfile="{dataset}/alignments/classif_by_coverage.out.log",
    benchmark:
        "{dataset}/alignments/classif_by_coverage.benchmark"
    conda:
        config.classif_by_coverage["conda"]
    threads: 1
    resources:
        disk_mb=1250,
        mem_mb=config.classif_by_coverage["mem"],
        runtime=config.classif_by_coverage["time"],
    params:
        THRESHOLD=config["classif_by_coverage"]["threshold"],
        GENE=config["classif_by_coverage"]["gene"],
        MIN_FRACT=config["classif_by_coverage"]["min_fraction_gene_covered"],
        MIN_COVERAGE=config["classif_by_coverage"]["min_avg_gene_coverage_depth"],
        CLASSIF_BY_COVERAGE=config.applications["classif_by_coverage"],
    shell:
        """
        {params.CLASSIF_BY_COVERAGE} --bed {input.CLASSIF_BED} --threshold {params.THRESHOLD} --gene "{params.GENE}" --min_fraction_gene_covered {params.MIN_FRACT} --min_avg_gene_coverage_depth {params.MIN_COVERAGE} --output {output.CLASSIF_CSV} -- {input.COVERAGE}    \
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
    log:
        outfile=cohortdir("coverage.out.log"),
        errfile=cohortdir("coverage.out.log"),
    benchmark:
        cohortdir("minority_variants.benchmark")
    conda:
        config.coverage["conda"]
    threads: config.coverage["threads"]
    resources:
        disk_mb=1250,
        mem_mb=config.coverage["mem"],
        runtime=config.coverage["time"],
    params:
        GATHER_COVERAGE=config.applications["gather_coverage"],
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
    log:
        outfile=cohortdir("minority_variants.out.log"),
        errfile=cohortdir("minority_variants.out.log"),
    benchmark:
        cohortdir("minority_variants.benchmark")
    conda:
        config.minor_variants["conda"]
    threads: config.minor_variants["threads"]
    resources:
        disk_mb=1250,
        mem_mb=config.minor_variants["mem"],
        runtime=config.minor_variants["time"],
    params:
        OUTDIR=cohortdir(""),
        NAMES=IDs,
        MIN_COVERAGE=config.minor_variants["min_coverage"],
        FREQUENCIES="--freqs" if config.minor_variants["frequencies"] else "",
        MINORITY_CALLER=config.applications["minority_freq"],
    shell:
        """
        {params.MINORITY_CALLER} -r {input.REF} -c {params.MIN_COVERAGE} -N {params.NAMES} -t {threads} -o {params.OUTDIR} {params.FREQUENCIES} {input.BAM} > >(tee {log.outfile}) 2> >(tee {log.errfile} >&2)
        """
