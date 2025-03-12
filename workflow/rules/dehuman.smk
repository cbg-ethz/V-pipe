from functools import partial


localrules:
    download_host_ref,


rule dh_reuse_alignreject:
    # "rely on aligner's output".
    # this rule re-use the rejected reads in align.smk (e.g. ngshmmalign's /alignments/rejects.sam)
    # (useful when in parallel with the main processing)
    input:
        reject_aln=(
            rules.hmm_align.output.reject_aln
            if config["general"]["aligner"] == "ngshmmalign"
            else temp_prefix("{dataset}/alignments/tmp_aln.sam")
        ),
    output:
        reject_1=temp_with_prefix("{dataset}/alignments/reject_R1.fastq.gz"),
        reject_2=temp_with_prefix("{dataset}/alignments/reject_R2.fastq.gz"),
    params:
        SAMTOOLS=config.applications["samtools"],
    log:
        outfile="{dataset}/alignments/reject.out.log",
        errfile="{dataset}/alignments/reject.err.log",
    conda:
        config.dehuman["conda"]
    benchmark:
        "{dataset}/alignments/reject.benchmark"
    group:
        "align"
    resources:
        disk_mb=1250,
        mem_mb=config.bwa_align["mem"],
        runtime=config.bwa_align["time"],
    threads: config.bwa_align["threads"]
    shell:
        """
        echo "Keep reject  -----------------------------------------------------"
        echo

        {params.SAMTOOLS} bam2fq -@ {threads} \
                                 -F 2 \
                                 -1 {output.reject_1} \
                                 -2 {output.reject_2} \
                                 {input.reject_aln} \
                                 2> >(tee {log.errfile} >&2)
        echo
        """


rule dh_redo_alignreject:
    # "redo rejects"
    # this rule re-does an alignment to get the reject
    # (useful if aligner's (temp) rejects have been deleted by now)
    input:
        global_ref=reference_file,
        ref_index=multiext(reference_file, *bwa_idx_ext),
        fastq=input_align_gz,
    output:
        tmp_aln=temp_with_prefix("{dataset}/alignments/dh_aln.sam"),
        reject_1=temp_with_prefix("{dataset}/alignments/reject_R1.fastq.gz"),
        reject_2=temp_with_prefix("{dataset}/alignments/reject_R2.fastq.gz"),
    params:
        BWA=config.applications["bwa"],
        SAMTOOLS=config.applications["samtools"],
    log:
        outfile="{dataset}/alignments/reject.out.log",
        errfile="{dataset}/alignments/reject.err.log",
    conda:
        config.dehuman["conda"]
    benchmark:
        "{dataset}/alignments/reject.benchmark"
    group:
        "dehuman"
    resources:
        disk_mb=1250,
        mem_mb=config.bwa_align["mem"],
        runtime=config.bwa_align["time"],
    threads: config.bwa_align["threads"]
    shell:
        """
        echo "Filter out virus' reads  -----------------------------------------"
        echo

        {params.BWA} mem -t {threads} \
                         -o {output.tmp_aln} \
                         {input.global_ref} {input.fastq} \
                         > {log.outfile} 2> >(tee {log.errfile} >&2)

        echo
        echo "Keep reject  -----------------------------------------------------"
        echo

        {params.SAMTOOLS} bam2fq -@ {threads} \
                                 -F 2 \
                                 -1 {output.reject_1} \
                                 -2 {output.reject_2} \
                                 {output.tmp_aln} \
                                 2> >(tee -a {log.errfile} >&2)
        echo
        """


rule download_host_ref:
    output:
        host_ref=config.dehuman["ref_host"],
    params:
        host_ref_url=config.dehuman["ref_host_url"],
    shell:
        """
        curl --output "{output.host_ref}" "{params.host_ref_url}"
        """


rule dh_hostalign:
    input:
        host_ref=config.dehuman["ref_host"],
        ref_index=multiext(config.dehuman["ref_host"], *bwa_idx_ext),
        reject_1=(
            rules.dh_redo_alignreject.output.reject_1
            if config["dehuman"]["catchup"]
            else rules.dh_reuse_alignreject.output.reject_1
        ),
        reject_2=(
            rules.dh_redo_alignreject.output.reject_2
            if config["dehuman"]["catchup"]
            else rules.dh_reuse_alignreject.output.reject_2
        ),
    output:
        host_aln=temp_with_prefix("{dataset}/alignments/host_aln.sam"),
    params:
        BWA=config.applications["bwa"],
    log:
        outfile="{dataset}/alignments/host_aln.out.log",
        errfile="{dataset}/alignments/host_aln.err.log",
    conda:
        config.dehuman["conda"]
    benchmark:
        "{dataset}/alignments/host_aln.benchmark"
    group:
        "dehuman"
    resources:
        disk_mb=1250,
        mem_mb=config.dehuman["mem"],
        runtime=config.dehuman["time"],
    threads: config.dehuman["threads"]
    shell:
        # create index if not exists:
        # test -f {input.ref_index} || {params.BWA} index {input.host_ref}
        """
        echo "Checking rejects against host's genome  --------------------------"

        {params.BWA} mem -t {threads} \
                         -o {output.host_aln}\
                         {input.host_ref} {input.reject_1} {input.reject_2} \
                         > {log.outfile} 2> >(tee {log.errfile} >&2)

        echo
        """


rule dh_filter:
    input:
        host_aln=temp_prefix("{dataset}/alignments/host_aln.sam"),
        R1=(
            partial(raw_data_file, pair=1)
            if config["dehuman"]["catchup"]
            else temp_prefix("{dataset}/extracted_data/R1.fastq")
        ),
        R2=(
            partial(raw_data_file, pair=2)
            if config["dehuman"]["catchup"]
            else temp_prefix("{dataset}/extracted_data/R2.fastq")
        ),
    output:
        filter_count="{dataset}/alignments/dehuman.count",
        filter_list=temp_with_prefix("{dataset}/alignments/dehuman.filter"),
        # TODO shift to pipe
        filtered_1=temp_with_prefix("{dataset}/raw_uploads/filtered_1.fastq.gz"),
        filtered_2=temp_with_prefix("{dataset}/raw_uploads/filtered_2.fastq.gz"),
    params:
        SAMTOOLS=config.applications["samtools"],
        remove_reads_script=cachepath(
            "../scripts/remove_reads_list.pl", executable=True, localsource=True
        ),
        # TODO shoo out the cats
        keep_host=int(config.dehuman["keep_host"]),
        sort_tmp=temp_prefix("{dataset}/alignments/host_sort.tmp"),
        host_aln_cram="{dataset}/alignments/host_aln.cram",
        # set to 1 to trigger matches with human genome (used for testing):
        F=2,
    log:
        outfile="{dataset}/raw_uploads/dehuman_filter.out.log",
        errfile="{dataset}/raw_uploads/dehuman_filter.err.log",
    conda:
        config.dehuman["conda"]
    benchmark:
        "{dataset}/raw_uploads/dehuman_filter.benchmark"
    group:
        "dehuman"
    resources:
        disk_mb=1250,
        mem_mb=config.dehuman["mem"],
        runtime=config.dehuman["time"],
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

        count=$({params.SAMTOOLS} view -@ {threads} -c -f {params.F} -F 2304 {input.host_aln} | tee {output.filter_count} 2> >(tee {log.errfile} >&2) )

        if (( count > 0 )); then
            echo
            echo "-----------------------------------------------------------------"
            echo "Needs special care: ${{count}} potential human reads found"
            echo "-----------------------------------------------------------------"
            echo
            echo "Removing identified host reads from raw reads -------------------"
            echo

            # get list
            {params.SAMTOOLS} view -@ {threads} \
                                   -f {params.F} \
                                   {input.host_aln} \
                                   | cut -f 1 > {output.filter_list} 2> >(tee -a {log.errfile} >&2)

            unpack_rawreads {input.R1:q} \
                   | {params.remove_reads_script} {output.filter_list} \
                   | gzip \
                   > {output.filtered_1} 2> >(tee -a {log.errfile} >&2) &

            unpack_rawreads {input.R2:q} \
                   | {params.remove_reads_script} {output.filter_list} \
                   | gzip \
                   > {output.filtered_2} 2> >(tee -a {log.errfile} >&2) &

            wait

            if (( {params.keep_host} )); then
                # keep the rejects for further analysis

                echo
                echo "Keeping host-aligned virus' rejects ------------------------------"
                echo


                # (we compress reference-less, because the reference size is larger
                # than the contaminant reads)
                rm -f '{params.sort_tmp}'.[0-9]*.bam
                FMT=cram,no_ref,use_bzip2,use_lzma,level=9,seqs_per_slice=1000000
                {params.SAMTOOLS} view -@ {threads} \
                                      -h -f 2 \
                                      {input.host_aln} \
                    | {params.SAMTOOLS} sort -@ {threads} \
                                    -T {params.sort_tmp} \
                                    --output-fmt ${{FMT}} \
                                    -o {params.host_aln_cram} \
                                    2> >(tee -a {log.errfile} >&2)

                echo
                echo "Compressing host-depleted raw reads ------------------------------"
                echo

                {params.SAMTOOLS} index -@ {threads} {params.host_aln_cram} 2> >(tee -a {log.errfile} >&2)
            fi
        else
            echo
            echo "No potential human reads found -----------------------------------"
            echo "Copy raw reads file"
            echo
            unpack_rawreads {input.R1} | gzip > {output.filtered_1} 2> >(tee -a {log.errfile} >&2) &
            unpack_rawreads {input.R2} | gzip > {output.filtered_2} 2> >(tee -a {log.errfile} >&2) &
            wait
            touch {output.filter_list}
        fi
        echo
        """


# TODO move cram compression into align.smk


rule dehuman:
    input:
        global_ref=reference_file,
        ref_index=multiext(reference_file, *bwa_idx_ext),
        filtered_1=rules.dh_filter.output.filtered_1,  # =temp_prefix("{dataset}/raw_uploads/filtered_1.fastq.gz"),
        filtered_2=rules.dh_filter.output.filtered_2,  # =temp_prefix("{dataset}/raw_uploads/filtered_2.fastq.gz"),
    output:
        cram_sam=temp_with_prefix("{dataset}/raw_uploads/dehuman.sam"),
        final_cram="{dataset}/raw_uploads/dehuman.cram",
        checksum="{dataset}/raw_uploads/dehuman.cram.%s" % config.general["checksum"],
    params:
        BWA=config.applications["bwa"],
        SAMTOOLS=config.applications["samtools"],
        checksum_type=config.general["checksum"],
        sort_tmp=temp_prefix("{dataset}/raw_uploads/dehuman.tmp"),
        # as a param to escape backslashes
        REGEXP=r"s{(?<=\t)([[:digit:]]:[[:upper:]]:[[:digit:]]:([ATCGN]+(\+[ATCGN]+)?|[[:digit:]]+))$}{BC:Z:\1}",
    log:
        outfile="{dataset}/raw_uploads/dehuman.out.log",
        errfile="{dataset}/raw_uploads/dehuman.err.log",
    conda:
        config.dehuman["conda"]
    benchmark:
        "{dataset}/raw_uploads/dehuman.benchmark"
    group:
        "dehuman"
    resources:
        disk_mb=1250,
        mem_mb=config.dehuman["mem"],
        runtime=config.dehuman["time"],
    threads: config.dehuman["threads"]
    shell:
        """
        echo "Compress filtered sequences --------------------------------------"
        echo

        {params.BWA} mem -t {threads} \
                         -C \
                         -o {output.cram_sam} \
                         {input.global_ref} {input.filtered_1} {input.filtered_2} \
                         > {log.outfile} 2> >(tee {log.errfile} >&2)

        # HACK handle incompatibilities between:
        #  - Illumina's 'bcl2fastq', which write arbitrary strings
        #  - 'bwa mem' which keep comments in the SAM file verbatim as in the FASTQ file
        #  - 'samtools' which expects comment to be properly marked as 'BC:Z:'
        #    as per SAM format specs
        REGEXP=\'{params.REGEXP}\'
        FMT=cram,embed_ref,use_bzip2,use_lzma,level=9,seqs_per_slice=1000000

        rm -f '{params.sort_tmp}'.[0-9]*.bam
        perl -p -e ${{REGEXP}} {output.cram_sam} \
              | {params.SAMTOOLS} sort -@ {threads} \
                                       -T {params.sort_tmp} \
                                       -M \
                                       --reference {input.global_ref} \
                                       --output-fmt ${{FMT}} \
                                       -o {output.final_cram} \
                                       2> >(tee -a {log.errfile} >&2)

        {params.checksum_type}sum {output.final_cram} > {output.checksum} 2> >(tee -a {log.errfile} >&2)

        echo
        echo DONE -------------------------------------------------------------
        echo
        """
