# TODO: how to set date?


configfile: "vpipeFetch.config.yaml"


rule all:
    input:
        "samples.tsv",


rule download_fastq:
    output:
        fname_info=temp("samples/{accession}/19700101/raw_data/info.tsv"),
    params:
        restart_times=3,
    log:
        outfile="logs/download.{accession}.out.log",
        errfile="logs/download.{accession}.err.log",
    conda:
        "envs/sratools.yaml"
    resources:
        mem_mb=5_000,
    threads: 6
    script:
        "scripts/download.py"


rule create_samples_tsv:
    input:
        expand(
            "samples/{accession}/19700101/raw_data/info.tsv",
            accession=config["samples"],
        ),
    output:
        "samples.tsv",
    shell:
        """
        cat {input} > {output}
        """
