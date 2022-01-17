from functools import partial
from glob import glob


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


def raw_data_file(wildcards, pair):
    for p in os.listdir('{dataset}/raw_data'.format(dataset=wildcards.dataset)):
        if re.search(r".*R{pair}\.(fastq\.gz|fastq|fq|fq\.gz)$".format(pair=pair), p):
            return os.path.join(wildcards.dataset, "raw_data", p)


def temp_with_prefix(p):
    return os.path.join(config.general["temp_prefix"], p)

        

rule dehuman:
    input:
        global_ref=reference_file,
        R1=partial(raw_data_file, pair=1),
        R2=partial(raw_data_file, pair=2),
    output:
        final_cram="{dataset}/raw_data/dehuman.cram",
        ref_aln=temp_with_prefix("{dataset}/ref_aln.sam"),
        filter=temp_with_prefix("{dataset}/dehuman.filter"),
        filtered_1=temp_with_prefix("{dataset}/raw_data/filtered_1.fasta.gz"),
        filtered_2=temp_with_prefix("{dataset}/raw_data/filtered_2.fasta.gz"),
        reject_1=temp_with_prefix("{dataset}/reject_R1.fastq.gz"),
        reject_2=temp_with_prefix("{dataset}/reject_R2.fastq.gz"),
        h38_aln=temp_with_prefix("{dataset}/h38_aln.sam"),
        h38_aln_cram=temp_with_prefix("{dataset}/h38_aln.cram"),
        cram_sam=temp_with_prefix("{dataset}/raw_data/cram.sam"),
    threads: config.dehuman["threads"]
    params:
        BWA=config.applications["bwa"],
        SAMTOOLS=config.applications["samtools"],
        HUMAN_GENOME=config.dehuman["ref_human"],
        remove_reads_script=cachepath(
            "../scripts/remove_reads_list.pl", executable=True, localsource=True
        ),
        stats="{dataset}/raw_data/dehuman.stats",
        # set to 1 to trigger matches with human genome (used for testing):
        F=2
    conda:
        config.dehuman["conda"]
    shell:
        """
        # create index if not exists:
        test -f {params.HUMAN_GENOME}.bwt || {params.BWA} index {params.HUMAN_GENOME}

        echo Filter out SARS-CoV-2 reads  -------------------------------------
        echo

        {params.BWA} mem -t {threads} \
                         -o {output.ref_aln} \
                         {input.global_ref} {input.R1} {input.R2}

        echo
        echo Keep reject  -----------------------------------------------------
        echo

        {params.SAMTOOLS} bam2fq -@ {threads} \
                                 -F 2 \
                                 -1 {output.reject_1} \
                                 -2 {output.reject_2} \
                                 {output.ref_aln}

        echo
        echo Checking rejects against Homo Sapiens ---------------------------

        {params.BWA} mem -t {threads} \
                         -o {output.h38_aln}\
                         {params.HUMAN_GENOME} {output.reject_1} {output.reject_2}

        echo
        echo Count aligned reads ---------------------------------------------
        echo

        count=$({params.SAMTOOLS} view -@ {threads} -c -f {params.F} {output.h38_aln})

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
                                   {output.h38_aln} \
                                   | cut -f 1 > {output.filter}

            # using zcat FILENAME.gz causes issues on Mac, see
            # https://serverfault.com/questions/570024/
            # redirection fixes this:
            zcat < {wildcards.dataset}/raw_data/*_R1.fastq.gz \
                   | {params.remove_reads_script} {output.filter} \
                   | gzip \
                   > {output.filtered_1} &

            zcat < {wildcards.dataset}/raw_data/*_R2.fastq.gz \
                   | {params.remove_reads_script} {output.filter} \
                   | gzip \
                   > {output.filtered_2} &

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
                                   -o {output.h38_aln_cram} \
                                   {output.h38_aln}

            echo
            echo Compressing human-depleted raw reads -----------------------------
            echo

            {params.SAMTOOLS} index -@ {threads} {output.h38_aln_cram}
            {params.SAMTOOLS} stats -@ {threads} {output.h38_aln_cram} > {params.stats}

        else
            echo
            echo No potential human reads found -----------------------------------
            echo Copy raw reads file
            echo
            cp {wildcards.dataset}/raw_data/*_R1.fastq.gz {output.filtered_1} &
            cp {wildcards.dataset}/raw_data/*_R2.fastq.gz {output.filtered_2} &
            wait
            touch {output.filter} 
            touch {output.h38_aln_cram}
        fi

        echo
        echo Compress filtered sequences --------------------------------------
        echo

        {params.BWA} mem -t {threads} \
                         -C \
                         -o {output.cram_sam} \
                         {input.global_ref} {output.filtered_1} {output.filtered_2}

        REGEXP=\'s{{(?<=\\t)([[:digit:]]:[[:upper:]]:[[:digit:]]:[ATCGN]+(\+[ATCGN]+))$}}{{BC:Z:\\1}}\'
        FMT=cram,embed_ref,use_bzip2,use_lzma,level=9,seqs_per_slice=1000000

        perl -p -e ${{REGEXP}} {output.cram_sam} \
              | {params.SAMTOOLS} sort -@ {threads} \
                                       -M \
                                       --reference {input.global_ref} \
                                       --output-fmt ${{FMT}} \
                                       -o {output.final_cram}
        echo
        echo DONE -------------------------------------------------------------
        echo
        """
