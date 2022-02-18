if config.upload["conda"]:

    localrules:
        prepare_upload,


def bcft_suffix():
    return (
        ""
        if config.upload["consensus"] == "majority"
        else f'_{config.upload["consensus"]}'
    )


rule prepare_upload:
    input:
        R1=partial(raw_data_file, pair=1),
        R2=partial(raw_data_file, pair=2),
        dehuman=[
            "{dataset}/raw_uploads/dehuman.cram",
            "{dataset}/raw_uploads/dehuman.cram.%s" % config.general["checksum"],
        ]
        if config.output["dehumanized_raw_reads"]
        else [],
        consensus_indels="{dataset}/references/consensus%s.bcftools.fasta"
        % bcft_suffix(),
        consensus_indels_chain="{dataset}/references/consensus%s.bcftools.chain"
        % bcft_suffix(),
        consensus_aligned="{dataset}/references/ref_majority_dels.fasta",  # % config.upload["consensus"],
        csum=[
            "{dataset}/references/consensus%(suf)s.bcftools.fasta.%(csum)s"
            % {"suf": bcft_suffix(), "csum": config.general["checksum"]},
            "{dataset}/references/ref_majority_dels.fasta.%(csum)s"
            % {"suf": bcft_suffix(), "csum": config.general["checksum"]},
        ]
        if config.upload["checksum"]
        else [],
    output:
        upload_prepared_touch="{dataset}/upload_prepared.touch",
    params:
        sample_id=ID,
        script=cachepath(config.upload["script"], executable=True),
        options=config.upload["options"],
    conda:
        # NOTE realpath is a gnu coreutils executable and not available out of the box. We need a conda environment anyway
        config.upload["conda"]
    resources:
        disk_mb=1000,
        mem_mb=config.upload["mem"],
        time_min=config.upload["time"],
    threads: config.upload["threads"]
    shell:
        """
        {params.script} {params.options} "{output.upload_prepared_touch}" "{params.sample_id}" "{wildcards.dataset}" {input}
        """


rule checksum:
    input:
        "{file}",
    output:
        "{file}.%s" % config.general["checksum"],
    params:
        checksum_type=config.general["checksum"],
    conda:
        config["upload"]["conda"]
    resources:
        disk_mb=10,
        mem_mb=config.checksum["mem"],
        time_min=config.checksum["time"],
    threads: 1
    shell:
        """
        {params.checksum_type}sum {input} > {output}
        """


ruleorder: dehuman > checksum
