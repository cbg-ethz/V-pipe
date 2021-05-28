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
        temp("{file}.{ext,(fastq|fq)}"),
    params:
        scratch="10000",
        mem=config.gunzip["mem"],
        time=config.gunzip["time"],
        GUNZIP=config.applications["gunzip"],
    log:
        outfile=temp("{file}_{ext}_gunzip.out.log"),
        errfile=temp("{file}_{ext}_gunzip.err.log"),
    threads: 1
    shell:
        """
        {params.GUNZIP} -c {input} > {output}
        """


rule extract:
    input:
        construct_input_fastq,
    output:
        temp("{dataset}/extracted_data/R{pair}.fastq"),
    params:
        scratch="2000",
        mem=config.extract["mem"],
        time=config.extract["time"],
    log:
        outfile="{dataset}/extracted_data/extract_R{pair}.out.log",
        errfile="{dataset}/extracted_data/extract_R{pair}.err.log",
    benchmark:
        "{dataset}/extracted_data/extract_R{pair}.benchmark"
    threads: 1
    shell:
        """
        cat {input} | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > {output} 2> >(tee {log.errfile} >&2)
        """


# 2. clipping
def len_cutoff(wildcards):
    parts = wildcards.dataset.split("/")
    patient_ID = parts[1]
    date = parts[2]
    patient_tuple = patient_record(patient_id=patient_ID, date=date)
    read_len = patient_dict[patient_tuple]
    len_cutoff = int(config.input["trim_percent_cutoff"] * read_len)
    return len_cutoff


if config.input["paired"]:

    rule preprocessing:
        input:
            R1="{dataset}/extracted_data/R1.fastq",
            R2="{dataset}/extracted_data/R2.fastq",
        output:
            R1gz="{dataset}/preprocessed_data/R1.fastq.gz",
            R2gz="{dataset}/preprocessed_data/R2.fastq.gz",
        params:
            scratch="2000",
            mem=config.preprocessing["mem"],
            time=config.preprocessing["time"],
            EXTRA=config.preprocessing["extra"],
            LEN_CUTOFF=len_cutoff,
            PRINSEQ=config.applications["prinseq"],
        log:
            outfile="{dataset}/preprocessed_data/prinseq.out.log",
            errfile="{dataset}/preprocessed_data/prinseq.err.log",
        conda:
            config.preprocessing["conda"]
        benchmark:
            "{dataset}/preprocessed_data/prinseq.benchmark"
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
            R1="{dataset}/extracted_data/R1.fastq",
        output:
            R1gz="{dataset}/preprocessed_data/R1.fastq.gz",
        params:
            scratch="2000",
            mem=config.preprocessing["mem"],
            time=config.preprocessing["time"],
            EXTRA=config.preprocessing["extra"],
            LEN_CUTOFF=len_cutoff,
            PRINSEQ=config.applications["prinseq"],
        log:
            outfile="{dataset}/preprocessed_data/prinseq.out.log",
            errfile="{dataset}/preprocessed_data/prinseq.err.log",
        conda:
            config.preprocessing["conda"]
        benchmark:
            "{dataset}/preprocessed_data/prinseq.benchmark"
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
        "{dataset}/extracted_data/R{pair}.fastq",
    output:
        "{dataset}/extracted_data/R{pair}_fastqc.html",
    params:
        scratch="2000",
        mem=config.fastqc["mem"],
        time=config.fastqc["time"],
        NOGROUP="--nogroup" if config.fastqc["no_group"] else "",
        OUTDIR="{dataset}/extracted_data",
        FASTQC=config.applications["fastqc"],
    log:
        outfile="{dataset}/extracted_data/R{pair}_fastqc.out.log",
        errfile="{dataset}/extracted_data/R{pair}_fastqc.err.log",
    conda:
        config.fastqc["conda"]
    threads: config.fastqc["threads"]
    shell:
        """
        {params.FASTQC} -o {params.OUTDIR} -t {threads} {params.NOGROUP} {input} 2> >(tee {log.errfile} >&2)
        """
