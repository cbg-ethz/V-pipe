if config.upload["local"]:

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
        R1=partial(raw_data_file, pair=1) if config.upload["orig_fastq"] else [],
        R2=partial(raw_data_file, pair=2) if config.upload["orig_fastq"] else [],
        orig_cram=[
            "{dataset}/raw_uploads/raw_reads.cram",
            "{dataset}/raw_uploads/raw_reads.cram.%s" % config.general["checksum"],
        ]
        if config.upload["orig_cram"]
        else [],
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
        consensus_aligned="{dataset}/references/ref_%s_dels.fasta"
        % config.upload["consensus"],
        csum=[
            "{dataset}/references/consensus%(suf)s.bcftools.fasta.%(csum)s"
            % {"suf": bcft_suffix(), "csum": config.general["checksum"]},
            "{dataset}/references/ref_majority_dels.fasta.%(csum)s"
            % {"suf": bcft_suffix(), "csum": config.general["checksum"]},
        ]
        if config.upload["checksum"]
        else [],
        frameshift_deletions_check="{dataset}/references/frameshift_deletions_check.tsv"
        if config.output["QA"]
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
        {params.script} {params.options} "{output.upload_prepared_touch}" "{params.sample_id}" "{wildcards.dataset}" {input:q}
        """


rule unfiltered_cram:
    input:
        global_ref=reference_file,
        ref_index=multiext(reference_file, *bwa_idx_ext),
        # NOTE not relying on rule extract - if done as a catchup these might not exist
        R1=partial(raw_data_file, pair=1),
        R2=partial(raw_data_file, pair=2),
    output:
        cram_sam=temp_with_prefix("{dataset}/raw_uploads/raw_reads.sam"),
        final_cram="{dataset}/raw_uploads/raw_reads.cram",
        checksum="{dataset}/raw_uploads/raw_reads.cram.%s" % config.general["checksum"],
    params:
        BWA=config.applications["bwa"],
        SAMTOOLS=config.applications["samtools"],
        checksum_type=config.general["checksum"],
        sort_tmp=temp_prefix("{dataset}/raw_uploads/raw_reads.tmp"),
    conda:
        config.dehuman["conda"]
    resources:
        disk_mb=1250,
        mem_mb=config.bwa_align["mem"],
        time_min=config.bwa_align["time"],
    threads: config.bwa_align["threads"]
    shell:
        """
        # using zcat FILENAME.gz causes issues on Mac, see
        # https://serverfault.com/questions/570024/
        # redirection fixes this:
        unpack_rawreads() {{
            for fq in "${{@}}"; do
                zcat -f < "${{fq}}"
            done
        }}

        echo "Compress un-filtered sequences -----------------------------------"
        echo

        {params.BWA} mem -t {threads} \
                         -C \
                         -o {output.cram_sam} \
                         {input.global_ref} \
                         <(unpack_rawreads {input.R1:q}) \
                         <(unpack_rawreads {input.R2:q})

        # HACK handle incompatibilities between:
        #  - Illumina's 'bcl2fastq', which write arbitrary strings
        #  - 'bwa mem' which keep comments in the SAM file verbatim as in the FASTQ file
        #  - 'samtools' which expects comment to be properly marked as 'BC:Z:'
        #    as per SAM format specs
        REGEXP=\'s{{(?<=\\t)([[:digit:]]:[[:upper:]]:[[:digit:]]:([ATCGN]+(\+[ATCGN]+)?|[[:digit:]]+))$}}{{BC:Z:\\1}}\'
        FMT=cram,embed_ref,use_bzip2,use_lzma,level=9,seqs_per_slice=1000000

        rm -f '{params.sort_tmp}'.[0-9]*.bam
        perl -p -e ${{REGEXP}} {output.cram_sam} \
              | {params.SAMTOOLS} sort -@ {threads} \
                                       -T {params.sort_tmp} \
                                       -M \
                                       --reference {input.global_ref} \
                                       --output-fmt ${{FMT}} \
                                       -o {output.final_cram}

        {params.checksum_type}sum {output.final_cram} > {output.checksum}

        echo
        echo DONE -------------------------------------------------------------
        echo
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


ruleorder: unfiltered_cram > checksum
ruleorder: dehuman > checksum
