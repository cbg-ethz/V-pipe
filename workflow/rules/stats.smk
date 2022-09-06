__author__ = "Susana Posada-Cespedes"
__author__ = "David Seifert"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"


localrules:
    stats,
    aggregate_stats,


rule allstats:
    input:
        "stats/coverage_intervals.tsv",
        "stats/read_counts.tsv",


rule alignment_coverage:
    input:
        BAM=expand(alignment_wildcard, dataset=datasets),
        TSV="variants/coverage.tsv",
    output:
        "stats/coverage_intervals.tsv",
    params:
        COVERAGE=config.alignment_coverage["coverage"],
        NAMES=IDs,
        EXTRACT_COVERAGE_INTERVALS=config.applications["extract_coverage_intervals"],
        ARRAYBASED=config.general["tsvbased"],
    log:
        outfile="stats/alignment_coverage.out.log",
        errfile="stats/alignment_coverage.out.log",
    conda:
        config.alignment_coverage["conda"]
    benchmark:
        "stats/alignment_coverage.benchmark"
    resources:
        disk_mb=1250,
        mem_mb=config.alignment_coverage["mem"],
        time_min=config.alignment_coverage["time"],
    threads: 1
    shell:
        """
        {params.EXTRACT_COVERAGE_INTERVALS} -b {params.ARRAYBASED} -cf {input.TSV} -c {params.COVERAGE} --no-shorah -N {params.NAMES} -o {output} {input.BAM} > >(tee {log.outfile}) 2>&1
        """


rule stats:
    input:
        R1=construct_input_fastq,
        R1_QC="{dataset}/preprocessed_data/R1.fastq.gz",
        BAM=alignment_wildcard,
    output:
        temp("{dataset}/read_counts_{pair,1}.tsv"),
    params:
        ID=ID,
        R1_temp=lambda wildcards: f"{wildcards.dataset}/preprocessed_data/temp.fastq",
        # int(4) is a workaround for a snakefmt bug:
        FACTOR=2 if config.input["paired"] else int(4),
        SAMTOOLS=config.applications["samtools"],
        GUNZIP=config.applications["gunzip"],
    conda:
        config.stats["conda"]
    shell:
        """
        SAMPLE_ID="{params.ID}"

        # Number of input reads
        LINECOUNT=$( cat {input.R1} | wc -l )
        let "INPUT=LINECOUNT / {params.FACTOR}"

        # Number of reads after QC
        # For portability reason not using zcat
        {params.GUNZIP} -c {input.R1_QC} > {params.R1_temp}
        LINECOUNT=$( cat {params.R1_temp} | wc -l )
        let "READCOUNT=LINECOUNT / {params.FACTOR}"
        rm {params.R1_temp}

        # Number of aligned reads
        ALNCOUNT=$( {params.SAMTOOLS} view {input.BAM} | wc -l )

        echo -e "${{SAMPLE_ID}}\t${{INPUT}}\t${{READCOUNT}}\t${{ALNCOUNT}}" > {output}
        """


rule aggregate_stats:
    input:
        expand("{dataset}/read_counts_1.tsv", dataset=datasets),
    output:
        "stats/read_counts.tsv",
    shell:
        """
        cat {input} > stats/temp
        echo -e "ID\tInput\tQC\tAlignments" | cat - stats/temp > {output}
        rm stats/temp
        """
