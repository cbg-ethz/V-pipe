import os

__author__ = "Susana Posada-Cespedes"
__author__ = "David Seifert"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"


# 1. extract
rule gunzip:
    input:
        "{file}.{ext}.gz",
    output:
        pipe("{file}.{ext,(fastq|fq)}"),
    params:
        GUNZIP=config.applications["gunzip"],
    log:
        outfile=temp("{file}_{ext}_gunzip.out.log"),
        errfile=temp("{file}_{ext}_gunzip.err.log"),
    benchmark:
        "{file}_{ext}_gunzip.benchmark"
    group:
        "extract"
    resources:
        disk_mb=1000,
        mem_mb=config.gunzip["mem"],
        runtime=config.gunzip["time"],
    threads: 1
    shell:
        """
        {params.GUNZIP} -c {input} > {output}
        """


rule extract:
    input:
        # TODO replace with raw_data_file
        construct_input_fastq,
    output:
        temp_with_prefix("{dataset}/extracted_data/R{pair}.fastq"),
    log:
        outfile="{dataset}/extracted_data/extract_R{pair}.out.log",
        errfile="{dataset}/extracted_data/extract_R{pair}.err.log",
    benchmark:
        "{dataset}/extracted_data/extract_R{pair}.benchmark"
    group:
        "extract"
    resources:
        disk_mb=32768,  # for large files sort stores its temp data on disk
        mem_mb=config.extract["mem"],
        runtime=config.extract["time"],
    threads: 1
    shell:
        # TODO replace with better dedicated software
        """
        cat {input:q} | paste - - - - | LC_ALL=C sort -s -k1,1 -t " " | tr "\t" "\n" > {output} 2> >(tee {log.errfile} >&2)
        """


# 2. clipping
def len_cutoff(wildcards):
    s_rec = guess_sample(wildcards.dataset)
    read_len = sample_table[s_rec].len
    len_cutoff = int(config.input["trim_percent_cutoff"] * read_len)
    return len_cutoff


if not config.general["preprocessor"] or config.general["preprocessor"] == "skip":

    rule preprocessing_skip:
        input:
            temp_prefix("{dataset}/extracted_data/{pair}.fastq"),
        output:
            "{dataset}/preprocessed_data/{pair}.fastq.gz",
        log:
            outfile="{dataset}/preprocessed_data/skip-{pair}.out.log",
            errfile="{dataset}/preprocessed_data/skip-{pair}.err.log",
        benchmark:
            "{dataset}/preprocessed_data/skip-{pair}.benchmark"
        resources:
            disk_mb=20000,
            mem_mb=config.preprocessing["mem"],
            runtime=config.preprocessing["time"],
        threads: 1
        shell:
            """
            echo "Skipping preprocessing and compressing merged/sorted fastq files as-is" > {log.outfile}

            gzip -c {input} > {output} 2> >(tee {log.errfile} >&2)
            """

elif config.input["paired"]:

    rule preprocessing:
        input:
            R1=temp_prefix("{dataset}/extracted_data/R1.fastq"),
            R2=temp_prefix("{dataset}/extracted_data/R2.fastq"),
        output:
            R1gz="{dataset}/preprocessed_data/R1.fastq.gz",
            R2gz="{dataset}/preprocessed_data/R2.fastq.gz",
        params:
            EXTRA=config.preprocessing["extra"],
            LEN_CUTOFF=len_cutoff,
            PRINSEQ=config.applications["prinseq"],
        log:
            outfile="{dataset}/preprocessed_data/prinseq.out.log",
            errfile="{dataset}/preprocessed_data/prinseq.err.log",
        conda:
            config.preprocessing["conda"]
        shadow:
            "minimal"
        benchmark:
            "{dataset}/preprocessed_data/prinseq.benchmark"
        resources:
            disk_mb=20000,
            mem_mb=config.preprocessing["mem"],
            runtime=config.preprocessing["time"],
        threads: 1
        shell:
            """
            echo "The length cutoff is: {params.LEN_CUTOFF}" > {log.outfile}

            {params.PRINSEQ} -fastq {input.R1} -fastq2 {input.R2} {params.EXTRA} -out_format 3 -out_good {wildcards.dataset}/preprocessed_data/R -out_bad null -min_len {params.LEN_CUTOFF} -log {log.outfile} 2> >(tee {log.errfile} >&2)

            # make sure that the lock held prinseq has been effectively released and propagated
            # on some networked shares this could otherwise lead to confusion or corruption
            if [[ "$OSTYPE" =~ ^linux ]]; then
                echo "Waiting for unlocks" >&2
                for U in {wildcards.dataset}/preprocessed_data/R_{{1,2}}.fastq; do
                    flock -x -o ${{U}} -c "echo ${{U}} unlocked >&2"
                done
            fi

            mv {wildcards.dataset}/preprocessed_data/R{{_,}}1.fastq
            mv {wildcards.dataset}/preprocessed_data/R{{_,}}2.fastq
            rm -f {wildcards.dataset}/preprocessed_data/R_?_singletons.fastq

            gzip {wildcards.dataset}/preprocessed_data/R1.fastq
            gzip {wildcards.dataset}/preprocessed_data/R2.fastq
            """

else:

    rule preprocessing_se:
        input:
            R1=temp_prefix("{dataset}/extracted_data/R1.fastq"),
        output:
            R1gz="{dataset}/preprocessed_data/R1.fastq.gz",
        params:
            EXTRA=config.preprocessing["extra"],
            LEN_CUTOFF=len_cutoff,
            PRINSEQ=config.applications["prinseq"],
        log:
            outfile="{dataset}/preprocessed_data/prinseq.out.log",
            errfile="{dataset}/preprocessed_data/prinseq.err.log",
        conda:
            config.preprocessing["conda"]
        shadow:
            "minimal"
        benchmark:
            "{dataset}/preprocessed_data/prinseq.benchmark"
        resources:
            disk_mb=10000,
            mem_mb=config.preprocessing["mem"],
            runtime=config.preprocessing["time"],
        threads: 1
        shell:
            """
            echo "The length cutoff is: {params.LEN_CUTOFF}" > {log.outfile}

            {params.PRINSEQ} -fastq {input.R1} {params.EXTRA} -out_format 3 -out_good {wildcards.dataset}/preprocessed_data/R -out_bad null -min_len {params.LEN_CUTOFF} {params.EXTRA} -log {log.outfile} 2> >(tee {log.errfile} >&2)

            # make sure that the lock held prinseq has been effectively released and propagated
            # on some network shares this could otherwise lead to confusion or corruption
            if [[ "$OSTYPE" =~ ^linux ]]; then
                echo "Waiting for unlocks" >&2
                for U in {wildcards.dataset}/preprocessed_data/R_{{1,2}}.fastq; do
                    flock -x -o ${{U}} -c "echo ${{U}} unlocked >&2"
                done
            fi

            mv {wildcards.dataset}/preprocessed_data/R{{,1}}.fastq

            gzip {wildcards.dataset}/preprocessed_data/R1.fastq
            """


# 3. QC reports


rule fastqc:
    input:
        temp_prefix("{dataset}/extracted_data/R{pair}.fastq"),
    output:
        "{dataset}/extracted_data/R{pair}_fastqc.html",
    params:
        NOGROUP="--nogroup" if config.fastqc["no_group"] else "",
        OUTDIR="{dataset}/extracted_data",
        FASTQC=config.applications["fastqc"],
    log:
        outfile="{dataset}/extracted_data/R{pair}_fastqc.out.log",
        errfile="{dataset}/extracted_data/R{pair}_fastqc.err.log",
    conda:
        config.fastqc["conda"]
    resources:
        disk_mb=2000,
        mem_mb=config.fastqc["mem"],
        runtime=config.fastqc["time"],
    threads: config.fastqc["threads"]
    shell:
        """
        {params.FASTQC} -o {params.OUTDIR} -t {threads} {params.NOGROUP} {input} 2> >(tee {log.errfile} >&2)
        """
