rule write_summary_json:
    input:
        consensus = lambda wildcard: [f for f in all_files if f.endswith("ref_majority_dels.fasta")],
        dehumanized = lambda wildcard: [f for f in all_files if f.endswith("dehuman.cram")],
        all_files = all_files
    output:
        "summary.zip",
    conda:
        config.report_sequences["conda"]
    shell:
        # needed {CONDA_PREFIX} on my dev setup on Mac + pyenv:
        """
        ${{CONDA_PREFIX}}/bin/python {VPIPE_BASEDIR}/scripts/report_sequences.py {input.consensus}
        """

rule dehuman:
    input:
        global_ref=reference_file,
        # R1=expand("{dataset}/preprocessed_data/R1.fastq.gz", dataset=datasets),
        # R2=expand("{dataset}/preprocessed_data/R2.fastq.gz", dataset=datasets),
        R1="{dataset}/preprocessed_data/R1.fastq.gz",
        R2="{dataset}/preprocessed_data/R2.fastq.gz",
    output:
        final_cram="{dataset}/raw_data/dehuman.cram",
    threads: config.dehuman["threads"]
    params:
        BWA=config.applications["bwa"],
        SAMTOOLS=config.applications["samtools"],
        HUMAN_GENOME=config.dehuman["ref_human"],
        remove_reads_script=cachepath(
            "../scripts/remove_reads_list.pl", executable=True, localsource=True
        ),
        ref_aln="{dataset}/ref_aln.sam",
        filter="{dataset}/dehuman.filter",
        h38_aln="{dataset}/h38_aln.sam",
        h38_aln_cram="{dataset}/h38_aln.cram",
        stats="{dataset}/raw_data/dehuman.stats",
        # set to 1 to trigger matches with human genome (used for testing):
        F=2,
        reject_1="{dataset}/reject_R1.fastq.gz",
        reject_2="{dataset}/reject_R2.fastq.gz",
        filtered_1="{dataset}/raw_data/filtered_1.fasta.gz",
        filtered_2="{dataset}/raw_data/filtered_2.fasta.gz",
        cram_sam="{dataset}/raw_data/cram.sam",
    conda:
        config.dehuman["conda"]
    shell:
        """
        # create index if not exists:
        test -f {params.HUMAN_GENOME}.bwt || {params.BWA} index {params.HUMAN_GENOME}

        echo Filter out SARS-CoV-2 reads  -------------------------------------
        echo

        {params.BWA} mem -t {threads} \
                         -o {params.ref_aln} \
                         {input.global_ref} {input.R1} {input.R2}

        echo
        echo Keep reject  -----------------------------------------------------
        echo

        {params.SAMTOOLS} bam2fq -@ {threads} \
                                 -F 2 \
                                 -1 {params.reject_1} \
                                 -2 {params.reject_2} \
                                 {params.ref_aln}

        echo
        echo Checking rejects against Homo Sapiens ---------------------------


        {params.BWA} mem -t {threads} \
                         -o {params.h38_aln}\
                         {params.HUMAN_GENOME} {params.reject_1} {params.reject_2}

        echo
        echo Count aligned reads ---------------------------------------------
        echo

        count=$({params.SAMTOOLS} view -@ {threads} -c -f {params.F} {params.h38_aln})

        if (( count > 0 )); then
            echo
            echo -----------------------------------------------------------------
            echo "Needs special care: ${{count}} potential human reads found"
            echo -----------------------------------------------------------------
            echo
            echo Removing identified human reads from raw reads ------------------
            echo

            # get list
            {params.SAMTOOLS} view -@ {threads} \
                                   -f {params.F} \
                                   {params.h38_aln} \
                                   | cut -f 1 > {params.filter}

            # using zcat FILENAME.gz causes issues on Mac, see
            # https://serverfault.com/questions/570024/
            # redirection fixes this:
            zcat < {wildcards.dataset}/raw_data/*_R1.fastq.gz \
                   | {params.remove_reads_script} {params.filter} \
                   | gzip \
                   > {params.filtered_1} &

            zcat < {wildcards.dataset}/raw_data/*_R2.fastq.gz \
                   | {params.remove_reads_script} {params.filter} \
                   | gzip \
                   > {params.filtered_2} &

            wait

            # keep the rejects for further analysis

            echo
            echo Keeping Human-aligned SARS-CoV-2 rejects -------------------------
            echo

            # (we compress reference-less, because the reference size is larger
            # than the contaminant reads)
            FMT=cram,no_ref,use_bzip2,use_lzma,level=9,seqs_per_slice=1000000
            {params.SAMTOOLS} sort -@ {threads} \
                                   -M \
                                   --output-fmt ${{FMT}} \
                                   -o {params.h38_aln_cram} \
                                   {params.h38_aln}

            echo
            echo Compressing human-depleted raw reads -----------------------------
            echo

            {params.SAMTOOLS} index -@ {threads} {params.h38_aln_cram}
            {params.SAMTOOLS} stats -@ {threads} {params.h38_aln_cram} > {params.stats}

        else
            echo
            echo No potential human reads found -----------------------------------
            echo Copy raw reads file
            echo
            cp {wildcards.dataset}/raw_data/*_R1.fastq.gz {params.filtered_1} &
            cp {wildcards.dataset}/raw_data/*_R2.fastq.gz {params.filtered_2} &
            wait
        fi

        echo
        echo Compress filtered sequences --------------------------------------
        echo

        {params.BWA} mem -t {threads} \
                         -C \
                         -o {params.cram_sam} \
                         {input.global_ref} {params.filtered_1} {params.filtered_2}

        REGEXP=\'s{{(?<=\\t)([[:digit:]]:[[:upper:]]:[[:digit:]]:[ATCGN]+(\+[ATCGN]+))$}}{{BC:Z:\\1}}\'
        FMT=cram,embed_ref,use_bzip2,use_lzma,level=9,seqs_per_slice=1000000

        perl -p -e ${{REGEXP}} {params.cram_sam} \
              | {params.SAMTOOLS} sort -@ {threads} \
                                       -M \
                                       --reference {input.global_ref} \
                                       --output-fmt ${{FMT}} \
                                       -o {output.final_cram}
        echo
        echo DONE -------------------------------------------------------------
        echo
        """
