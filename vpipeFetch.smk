# TODO: how to set date?


configfile: "vpipeFetch.config.yaml"


rule all:
    input:
        "samples.tsv"


rule download:
    output:
        temp("samples/{accession}/19700101/raw_data/info.tsv")
    conda:
        "envs/sratools.yaml"
    log:
        "logs/download.{accession}.log"
    shell:
        """
        # TODO: make all of this more robust
        outdir="$(dirname {output[0]})"

        # download FastQ
        fasterq-dump \
            --outdir "$outdir" \
            {wildcards.accession} \
            &> "{log}"

        # make V-pipe recognize output files
        rename 's/_1/_R1/' "$outdir"/*.fastq
        rename 's/_2/_R2/' "$outdir"/*.fastq

        # get fastq filename
        fname=$(ls "$outdir"/*.fastq | head -1)
        echo "Selected filename: $fname" &>> "{log}"

        # get read length
        read_length=$(bioawk -c fastx '{{ bases += length($seq); count++ }} END{{print int(bases/count)}}' $fname)
        echo "Computed read length: $read_length" &>> "{log}"

        # store in info file
        echo -e "{wildcards.accession}\t19700101\t$read_length" > {output}
        """


rule create_samples_tsv:
    input:
        expand(
            "samples/{accession}/19700101/raw_data/info.tsv",
            accession=config["samples"]
        )
    output:
        "samples.tsv"
    shell:
        """
        cat {input} > {output}
        """
