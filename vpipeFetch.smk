# TODO: how to set date?


configfile: "vpipeFetch.config.yaml"


rule all:
    input:
        expand(
            "samples/{accession}/19700101/raw_data/info.csv",
            accession=config["samples"]
        )


rule download:
    output:
        "samples/{accession}/19700101/raw_data/info.csv"
    conda:
        "envs/sratools.yaml"
    log:
        "logs/download.{accession}.log"
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

        # mark finished download
        touch {output}
        """
