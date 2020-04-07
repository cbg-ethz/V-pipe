# TODO: how to set date?


configfile: "vpipeFetch.config.yaml"


rule all:
    input:
        expand(
            "samples/{accession}/19700101/raw_data/{accession}.fastq",
            accession=config["samples_se"]
        ),
        expand(
            "samples/{accession}/19700101/raw_data/{accession}_R1.fastq",
            accession=config["samples_pe"]
        )


rule download_se:
    output:
        "samples/{accession}/19700101/raw_data/{accession}.fastq",
    conda:
        "envs/sratools.yaml"
    log:
        "logs/fasterq-dump.{accession}.log"
    shell:
        """
        outdir="$(dirname {output[0]})"

        fasterq-dump \
            --outdir "$outdir" \
            {wildcards.accession} \
            &> "{log}"
        """


rule download_pe:
    output:
        "samples/{accession}/19700101/raw_data/{accession}_R1.fastq",
        "samples/{accession}/19700101/raw_data/{accession}_R2.fastq"
    conda:
        "envs/sratools.yaml"
    log:
        "logs/fasterq-dump.{accession}.log"
    shell:
        """
        outdir="$(dirname {output[0]})"

        fasterq-dump \
            --outdir "$outdir" \
            {wildcards.accession} \
            &> "{log}"

        # make V-pipe recognize output files
        # TODO: make this more robust
        rename 's/_1/_R1/' "$outdir"/*.fastq
        rename 's/_2/_R2/' "$outdir"/*.fastq
        """
