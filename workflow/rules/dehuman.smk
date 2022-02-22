from functools import partial


rule dh_reuse_alignreject:
    # "rely on aligner's output".
    # this rule re-use the rejected reads in align.smk (e.g. ngshmmalign's /alignments/rejects.sam)
    # (useful when in parallel with the main processing)
    input:
        reject_aln=rules.hmm_align.output.reject_aln
        if config["general"]["aligner"] == "ngshmmalign"
        else temp_prefix("{dataset}/alignments/tmp_aln.sam"),
    output:
        reject_1=temp_with_prefix("{dataset}/alignments/reject_R1.fastq.gz"),
        reject_2=temp_with_prefix("{dataset}/alignments/reject_R2.fastq.gz"),
    params:
        SAMTOOLS=config.applications["samtools"],
    conda:
        config.dehuman["conda"]
    group:
        "align"
    resources:
        disk_mb=1250,
        mem_mb=config.bwa_align["mem"],
        time_min=config.bwa_align["time"],
    threads: config.bwa_align["threads"]
    shell:
        """
        echo "Keep reject  -----------------------------------------------------"
        echo

        {params.SAMTOOLS} bam2fq -@ {threads} \
                                 -F 2 \
                                 -1 {output.reject_1} \
                                 -2 {output.reject_2} \
                                 {input.reject_aln}
        echo
        """


rule dh_redo_alignreject:
    # "redo rejects"
    # this rule re-does an alignment to get the reject
    # (useful if aligner's (temp) rejects have been deleted by now)
    input:
        global_ref=reference_file,
        ref_index="{}.bwt".format(reference_file),
        fastq=input_align_gz,
    output:
        tmp_aln=temp_with_prefix("{dataset}/alignments/dh_aln.sam"),
        reject_1=temp_with_prefix("{dataset}/alignments/reject_R1.fastq.gz"),
        reject_2=temp_with_prefix("{dataset}/alignments/reject_R2.fastq.gz"),
    params:
        BWA=config.applications["bwa"],
        SAMTOOLS=config.applications["samtools"],
    conda:
        config.dehuman["conda"]
    group:
        "dehuman"
    resources:
        disk_mb=1250,
        mem_mb=config.bwa_align["mem"],
        time_min=config.bwa_align["time"],
    threads: config.bwa_align["threads"]
    shell:
        """
        echo "Filter out virus' reads  -----------------------------------------"
        echo

        {params.BWA} mem -t {threads} \
                         -o {output.tmp_aln} \
                         {input.global_ref} {input.fastq}

        echo
        echo "Keep reject  -----------------------------------------------------"
        echo

        {params.SAMTOOLS} bam2fq -@ {threads} \
                                 -F 2 \
                                 -1 {output.reject_1} \
                                 -2 {output.reject_2} \
                                 {output.tmp_aln}
        echo
        """


rule dh_hostalign:
    input:
        host_ref=config.dehuman["ref_host"],
        ref_index="{}.bwt".format(config.dehuman["ref_host"]),
        reject_1=rules.dh_redo_alignreject.output.reject_1
        if config["dehuman"]["catchup"]
        else rules.dh_reuse_alignreject.output.reject_1,
        reject_2=rules.dh_redo_alignreject.output.reject_2
        if config["dehuman"]["catchup"]
        else rules.dh_reuse_alignreject.output.reject_2,
    output:
        host_aln=temp_with_prefix("{dataset}/alignments/host_aln.sam"),
    params:
        BWA=config.applications["bwa"],
    conda:
        config.dehuman["conda"]
    group:
        "dehuman"
    resources:
        disk_mb=1250,
        mem_mb=config.bwa_align["mem"],
        time_min=config.bwa_align["time"],
    threads: config.bwa_align["threads"]
    shell:
        # create index if not exists:
        # test -f {input.ref_index} || {params.BWA} index {input.host_ref}
        """
        echo "Checking rejects against Homo Sapiens ---------------------------"

        {params.BWA} mem -t {threads} \
                         -o {output.host_aln}\
                         {input.host_ref} {input.reject_1} {input.reject_2}

        echo
        """


rule dh_filter:
    input:
        host_aln=temp_prefix("{dataset}/alignments/host_aln.sam"),
        R1=partial(raw_data_file, pair=1),
        R2=partial(raw_data_file, pair=2),
    output:
        filter_count="{dataset}/alignments/dehuman.count",
        filter_list=temp_with_prefix("{dataset}/alignments/dehuman.filter"),
        # TODO shift to pipe
        filtered_1=temp_with_prefix("{dataset}/raw_uploads/filtered_1.fasta.gz"),
        filtered_2=temp_with_prefix("{dataset}/raw_uploads/filtered_2.fasta.gz"),
    params:
        SAMTOOLS=config.applications["samtools"],
        remove_reads_script=cachepath(
            "../scripts/remove_reads_list.pl", executable=True, localsource=True
        ),
        # TODO shoo out the cats
        keep_host=int(config.dehuman["keep_host"]),
        host_aln_cram="{dataset}/alignments/host_aln.cram",
        # set to 1 to trigger matches with human genome (used for testing):
        F=2,
    conda:
        config.dehuman["conda"]
    group:
        "dehuman"
    resources:
        disk_mb=1250,
        mem_mb=config.dehuman["mem"],
        time_min=config.dehuman["time"],
    threads: config.dehuman["threads"]
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

        echo
        echo "Count aligned reads ---------------------------------------------"
        echo

        count=$({params.SAMTOOLS} view -@ {threads} -c -f {params.F} {input.host_aln})
        echo "${{count}}" > {output.filter_count}

        if (( count > 0 )); then
            echo
            echo "-----------------------------------------------------------------"
            echo "Needs special care: ${{count}} potential human reads found"
            echo "-----------------------------------------------------------------"
            echo
            echo "Removing identified human reads from raw reads ------------------"
            echo

            # get list
            {params.SAMTOOLS} view -@ {threads} \
                                   -f {params.F} \
                                   {input.host_aln} \
                                   | cut -f 1 > {output.filter_list}

            unpack_rawreads {input.R1} \
                   | {params.remove_reads_script} {output.filter_list} \
                   | gzip \
                   > {output.filtered_1} &

            unpack_rawreads {input.R2} \
                   | {params.remove_reads_script} {output.filter_list} \
                   | gzip \
                   > {output.filtered_2} &

            wait

            if (( {params.keep_host} )); then
                # keep the rejects for further analysis

                echo
                echo "Keeping Human-aligned virus' rejects -----------------------------"
                echo

                # (we compress reference-less, because the reference size is larger
                # than the contaminant reads)
                FMT=cram,no_ref,use_bzip2,use_lzma,level=9,seqs_per_slice=1000000
                {params.SAMTOOLS} sort -@ {threads} \
                                    -M \
                                    --output-fmt ${{FMT}} \
                                    -o {params.host_aln_cram} \
                                    {input.host_aln}

                echo
                echo "Compressing human-depleted raw reads -----------------------------"
                echo

                {params.SAMTOOLS} index -@ {threads} {params.host_aln_cram}
            fi
        else
            echo
            echo "No potential human reads found -----------------------------------"
            echo "Copy raw reads file"
            echo
            # HACK cheating as currently we only use gziped data
            #unpack_rawreads {input.R1} | gzip > {output.filtered_1} &
            #unpack_rawreads {input.R2} | gzip > {output.filtered_2} &
            cat {input.R1} > {output.filtered_1} &
            cat {input.R2} > {output.filtered_2} &
            wait
            touch {output.filter_list}
        fi
        echo
        """


# TODO move cram compression into align.smk


rule dehuman:
    input:
        global_ref=reference_file,
        filtered_1=rules.dh_filter.output.filtered_1,  # =temp_prefix("{dataset}/raw_uploads/filtered_1.fasta.gz"),
        filtered_2=rules.dh_filter.output.filtered_2,  # =temp_prefix("{dataset}/raw_uploads/filtered_2.fasta.gz"),
    output:
        cram_sam=temp_with_prefix("{dataset}/raw_uploads/dehuman.sam"),
        final_cram="{dataset}/raw_uploads/dehuman.cram",
        checksum="{dataset}/raw_uploads/dehuman.cram.%s" % config.general["checksum"],
    params:
        BWA=config.applications["bwa"],
        SAMTOOLS=config.applications["samtools"],
        checksum_type=config.general["checksum"],
    conda:
        config.dehuman["conda"]
    group:
        "dehuman"
    resources:
        disk_mb=1250,
        mem_mb=config.bwa_align["mem"],
        time_min=config.bwa_align["time"],
    threads: config.bwa_align["threads"]
    shell:
        """
        echo "Compress filtered sequences --------------------------------------"
        echo

        {params.BWA} mem -t {threads} \
                         -C \
                         -o {output.cram_sam} \
                         {input.global_ref} {input.filtered_1} {input.filtered_2}

        # HACK handle incompatibilities between:
        #  - Illumina's 'bcl2fastq', which write arbitrary strings
        #  - 'bwa mem' which keep comments in the SAM file verbatim as in the FASTQ file
        #  - 'samtools' which expects comment to be properly marked as 'BC:Z:'
        #    as per SAM format specs
        REGEXP=\'s{{(?<=\\t)([[:digit:]]:[[:upper:]]:[[:digit:]]:([ATCGN]+(\+[ATCGN]+)?|[[:digit:]]+))$}}{{BC:Z:\\1}}\'
        FMT=cram,embed_ref,use_bzip2,use_lzma,level=9,seqs_per_slice=1000000

        perl -p -e ${{REGEXP}} {output.cram_sam} \
              | {params.SAMTOOLS} sort -@ {threads} \
                                       -M \
                                       --reference {input.global_ref} \
                                       --output-fmt ${{FMT}} \
                                       -o {output.final_cram}

        {params.checksum_type}sum {output.final_cram} > {output.checksum}

        echo
        echo DONE -------------------------------------------------------------
        echo
        """
