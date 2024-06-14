import os

__author__ = "Susana Posada-Cespedes"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"


def window_length1(wildcards):
    s_rec = guess_sample(wildcards.dataset)
    read_len = sample_table[s_rec].len
    aux = int((read_len * 4 / 5 + config.snv["shift"]) / config.snv["shift"])
    return str(aux * config.snv["shift"])


def shift1(wildcards):
    s_rec = guess_sample(wildcards.dataset)
    read_len = sample_table[s_rec].len
    aux = int((read_len * 4 / 5 + config.snv["shift"]) / config.snv["shift"])
    return str(aux)


# 1. Call single nucleotide variants
rule coverage_intervals:
    input:
        BAM=alignment_wildcard,
        TSV="{dataset}/alignments/coverage.tsv.gz",
    output:
        temp("{dataset}/variants/coverage_intervals.tsv"),
    params:
        NAME=ID,
        WINDOW_LEN=window_length1,
        SHIFT=shift1,
        COVERAGE=config.coverage_intervals["coverage"],
        OVERLAP="" if config.coverage_intervals["overlap"] else "-cf $TMPTSV",
        LIBERAL="-e" if config.coverage_intervals["liberal"] else "",
        EXTRACT_COVERAGE_INTERVALS=config.applications["extract_coverage_intervals"],
        GUNZIP=config.applications["gunzip"],
        ARRAYBASED=config.general["tsvbased"],
    log:
        outfile="{dataset}/variants/coverage_intervals.out.log",
        errfile="{dataset}/variants/coverage_intervals.out.log",
    conda:
        config.coverage_intervals["conda"]
    benchmark:
        "{dataset}/variants/coverage_intervals.benchmark"
    group:
        "snv"
    resources:
        disk_mb=1250,
        mem_mb=config.coverage_intervals["mem"],
        runtime=config.coverage_intervals["time"],
    threads: config.coverage_intervals["threads"]
    shell:
        """
        TMPTSV=$(mktemp -t XXXXXXXX_cov.tsv)
        TMPINT=$(mktemp -t XXXXXXXX_int.tsv)
        {params.GUNZIP} -c {input.TSV} | cut -f'2-' > $TMPTSV
        mkdir -p "$(dirname "{output}")"
        {params.EXTRACT_COVERAGE_INTERVALS} -c "{params.COVERAGE}" -w "{params.WINDOW_LEN}" -s "{params.SHIFT}" -N "{params.NAME}" {params.LIBERAL} -b {params.ARRAYBASED} {params.OVERLAP} -t "{threads}" -o $TMPINT "{input.BAM}" > >(tee {log.outfile}) 2>&1
        cat $TMPINT >> "{log.outfile}"
        read name intervals < $TMPINT
        IFS=',' read -r -a interarray <<< "$intervals"
        printf "%s\n" "${{interarray[@]}}" > "{output}"
        rm $TMPTSV $TMPINT
        """


def read_len(wildcards):
    s_rec = guess_sample(wildcards.dataset)
    read_len = sample_table[s_rec].len
    return read_len


rule snv:
    input:
        REF=(
            cohortdir("cohort_consensus.fasta")
            if config.snv["consensus"]
            else reference_file
        ),
        BAM=alignment_wildcard,
        TSV="{dataset}/variants/coverage_intervals.tsv",
    output:
        CSV="{dataset}/variants/SNVs/snvs.csv",
        VCF="{dataset}/variants/SNVs/snvs.vcf",
    params:
        READ_LEN=read_len,
        ALPHA=config.snv["alpha"],
        IGNORE_INDELS="--ignore_indels" if config.snv["ignore_indels"] else "",
        COVERAGE=config.snv["coverage"],
        SHIFT=config.snv["shift"],
        KEEP_FILES="true" if config.snv["keep_files"] else "false",
        seed="-R 42",
        WORK_DIR="{dataset}/variants/SNVs",
        LOCALSCRATCH=config.snv["localscratch"],
        SHORAH=config.applications["shorah"],
        POSTHRESH=config.snv["posterior_threshold"],
        COVINT=config.coverage_intervals["coverage"],
        BCFTOOLS=config.applications["bcftools"],
    log:
        outfile="{dataset}/variants/SNVs/shorah.out.log",
        errfile="{dataset}/variants/SNVs/shorah.err.log",
    conda:
        config.snv["conda"]
    benchmark:
        "{dataset}/variants/SNVs/shorah.benchmark"
    group:
        "snv"
    resources:
        disk_mb=1250,
        mem_mb=config.snv["mem"],
        runtime=config.snv["time"],
    threads: config.snv["threads"]
    shell:
        """
        let "WINDOW_SHIFTS=({params.READ_LEN} * 4/5 + {params.SHIFT}) / {params.SHIFT}"
        let "WINDOW_LEN=WINDOW_SHIFTS * {params.SHIFT}"

        echo "Windows are shifted by: ${{WINDOW_SHIFTS}} bp" > {log.outfile}
        echo "The window length is: ${{WINDOW_LEN}} bp" >> {log.outfile}

        # Get absolute path for input files
        CWD=${{PWD}}
        BAM=${{PWD}}/{input.BAM}
        REF={input.REF}; [[ ${{REF}} =~ ^/ ]] || REF=${{PWD}}/${{REF}}
        OUTFILE=${{PWD}}/{log.outfile}
        ERRFILE=${{PWD}}/{log.errfile}
        WORK_DIR=${{PWD}}/{params.WORK_DIR}

        if [[ -n "{params.LOCALSCRATCH}" ]]; then
            # put input files in localscratch
            rsync -aq "${{BAM}}" "${{REF}}" "{params.LOCALSCRATCH}/"
            BAM="{params.LOCALSCRATCH}/${{BAM##*/}}"
            REF="{params.LOCALSCRATCH}/${{REF##*/}}"
        fi

        # Run ShoRAH in each of the predetermined regions (regions with sufficient coverage)
        LINE_COUNTER=0
        FILES=( )
        FILES_VCF=( )
        while read -r region || [[ -n ${{region}} ]]
        do
            echo "Running ShoRAH on region: ${{region}}" >> $OUTFILE
            (( ++LINE_COUNTER ))
            # Create directory for running ShoRAH in a corresponding region (if doesn't exist)
            DIR=${{WORK_DIR}}/REGION_${{LINE_COUNTER}}
            if [[ ! -d "${{DIR}}" ]]; then
                echo "Creating directory ${{DIR}}" >> $OUTFILE
                mkdir -p ${{DIR}}
            else
                # Results from previous runs
                if [[ {params.KEEP_FILES} == "true" ]]; then
                    DIR_DST=${{WORK_DIR}}/old
                    echo "Moving results from a previous run to ${{DIR_DST}}" >> $OUTFILE
                    rm -rf ${{DIR_DST}}/REGION_${{LINE_COUNTER}}
                    mkdir -p ${{DIR_DST}}
                    mv -f ${{DIR}} ${{DIR_DST}}
                    mkdir -p ${{DIR}}
                fi
            fi
            # Change to the directory where ShoRAH is to be executed
            if [[ -z "{params.LOCALSCRATCH}" ]]; then
                cd ${{DIR}}
            else
                # special case: go in local scratch instead
                rsync -aq "${{DIR}}" "{params.LOCALSCRATCH}/"
                mkdir -p "{params.LOCALSCRATCH}/REGION_${{LINE_COUNTER}}"
                cd "{params.LOCALSCRATCH}/REGION_${{LINE_COUNTER}}"
            fi

            # NOTE: Execution command for ShoRAH2 valid from v1.99.0 and above
            {params.SHORAH} -t {threads} -a {params.ALPHA} -w ${{WINDOW_LEN}} -x 100000 {params.IGNORE_INDELS} -p {params.POSTHRESH} -c {params.COVERAGE} -r ${{region}} {params.seed} -b ${{BAM}} -f ${{REF}} >> $OUTFILE 2> >(tee -a $ERRFILE >&2)
            if [[ -n "{params.LOCALSCRATCH}" ]]; then
                # copyback from localscratch
                rsync -auq "{params.LOCALSCRATCH}/REGION_${{LINE_COUNTER}}" "${{WORK_DIR}}"
                cd ${{DIR}}
            fi
            if [[ -s ${{DIR}}/reads.fas && -f ${{DIR}}/snv/SNVs_0.010000_final.csv ]]; then
                # Non empty reads: run had enough data and should have produced SNVs
                {params.BCFTOOLS} view ${{DIR}}/snv/SNVs_0.010000_final.vcf -Oz -o ${{DIR}}/snv/SNVs_0.010000_final.vcf.gz
                {params.BCFTOOLS} index ${{DIR}}/snv/SNVs_0.010000_final.vcf.gz
                FILES+=("${{DIR}}/snv/SNVs_0.010000_final.csv")
                FILES_VCF+=("${{DIR}}/snv/SNVs_0.010000_final.vcf.gz")
            elif (( {params.COVINT} == 0 && LINE_COUNTER == 1 )) && [[ -f ${{DIR}}/reads.fas && ( ! -s ${{DIR}}/reads.fas ) ]]; then
                # if we have disabled coverage intervales entirely, the first and only line might have no reads
                # (e.g.: in negative controls )

                echo "No reads while coverage intervals disabled (possible negative control sample)" 2> >(tee -a $ERRFILE >&2)
                cd ${{CWD}}
                (( --LINE_COUNTER )) || true # Strict mode : (( 0 )) = fail
                break
            else
                echo "ERROR: unsuccesful execution of ShoRAH" 2> >(tee -a $ERRFILE >&2)
                exit 1
            fi

            # Change back to working directory
            cd ${{CWD}}
        done < {input.TSV}

        # Aggregate results from different regions
        if (( ${{#FILES[@]}} )); then
            echo "Intermediate csv files: ${{FILES[*]}}" >> {log.outfile}
            echo "Intermediate vcf files: ${{FILES_VCF[*]}}" >> {log.outfile}
            (head -n 1 "${{FILES[0]}}"; tail -q -n +2 "${{FILES[@]}}" | sort -t, -nk2) > {output.CSV}
            {params.BCFTOOLS} concat -o ${{WORK_DIR}}/snvs_tmp.vcf "${{FILES_VCF[@]}}"
            {params.BCFTOOLS} sort ${{WORK_DIR}}/snvs_tmp.vcf  -o {output.VCF}
            rm -f ${{WORK_DIR}}/snvs_tmp.vcf
        elif (( LINE_COUNTER )); then
            echo "ERROR: unsuccesful execution of ShoRAH" 2> >(tee -a {log.errfile} >&2)
            exit 1
        else
            echo "No alignment region reports sufficient coverage" >> {log.outfile}
            touch {output.CSV}
            touch {output.VCF}
        fi
        """


rule samtools_index:
    input:
        "{file}.fasta",
    output:
        "{file}.fasta.fai",
    params:
        SAMTOOLS=config.applications["samtools"],
    log:
        outfile="{file}_samtools_index.out.log",
        errfile="{file}_samtools_index.err.log",
    conda:
        config.samtools_index["conda"]
    resources:
        disk_mb=2000,
        mem_mb=config.samtools_index["mem"],
        runtime=config.samtools_index["time"],
    benchmark:
        "{file}_samtools_index.benchmark"
    shell:
        """
        {params.SAMTOOLS} faidx {input} -o {output} > {log.outfile} 2> >(tee -a {log.errfile} >&2)
        """


rule lofreq:
    input:
        REF=(
            cohortdir("cohort_consensus.fasta")
            if config.lofreq["consensus"]
            else reference_file
        ),
        REF_IDX=(
            cohortdir("cohort_consensus.fasta.fai")
            if config.lofreq["consensus"]
            else "%s.fai" % reference_file
        ),
        BAM=alignment_wildcard,
    output:
        BAM="{dataset}/variants/SNVs/REF_aln_indelqual.bam",
        SNVs="{dataset}/variants/SNVs/snvs.vcf",
    params:
        OUTDIR="{dataset}/variants/SNVs",
        EXTRA=config.lofreq["extra"],
        SAMTOOLS=config.applications["samtools"],
        LOFREQ=config.applications["lofreq"],
    log:
        outfile="{dataset}/variants/SNVs/lofreq.out.log",
        errfile="{dataset}/variants/SNVs/lofreq.err.log",
    conda:
        config.lofreq["conda"]
    benchmark:
        "{dataset}/variants/SNVs/lofreq.benchmark"
    resources:
        disk_mb=2000,
        mem_mb=config.lofreq["mem"],
        runtime=config.lofreq["time"],
    shell:
        """
        # Add qualities to indels
        {params.LOFREQ} indelqual --dindel -f {input.REF} -o {output.BAM} --verbose {input.BAM} > {log.outfile} 2> >(tee {log.errfile} >&2)
        # Index bam file
        {params.SAMTOOLS} index {output.BAM} 2> >(tee {log.errfile} >&2)

        # Run Lofreq
        echo "Running LoFreq" >> {log.outfile}
        {params.LOFREQ} call {params.EXTRA} --call-indels -f {input.REF} -o {output.SNVs} --verbose {output.BAM} >> {log.outfile} 2> >(tee -a {log.errfile} >&2)
        """


rule paired_end_read_merger:
    input:
        fname_bam=alignment_wildcard,
        fname_ref=(
            cohortdir("cohort_consensus.fasta")
            if config.viloca["consensus"]
            else reference_file
        ),
        fname_ref_idx=(
            cohortdir("cohort_consensus.fasta")
            if config.viloca["consensus"]
            else reference_file
        )
        + ".fai",
    output:
        fname_bam_merged="{dataset}/alignment/REF_aln.merged.bam",
        fname_sam_merged=temp_with_prefix("{dataset}/alignment/REF_aln.merged.sam"),
        fname_sam=temp_with_prefix("{dataset}/alignment/REF_aln.sam"),
        fname_sam_nonmerged="{dataset}/alignment/REF_aln.nonmerged.sam",
        fname_sam_sort=temp_with_prefix("{dataset}/alignment/REF_aln.sort.sam"),
    params:
        SAMTOOLS=config.applications["samtools"],
        PAIRED_END_READ_MERGER=config.applications["paired_end_read_merger"],
        sort_tmp=temp_prefix("{dataset}.tmp"),
    log:
        outfile="{dataset}/alignment/paired_end_read_merger.out.log",
        errfile="{dataset}/alignment/paired_end_read_merger.err.log",
    conda:
        config.paired_end_read_merger["conda"]
    shell:
        """
        ## Preparation
        {params.SAMTOOLS} view -h -T {input.fname_ref} -t {input.fname_ref_idx} {input.fname_bam} -o {output.fname_sam} > {log.outfile} 2> >(tee {log.errfile} >&2)
        ## sort accrording to QNAME
        rm -f '{params.sort_tmp}'.[0-9]*.bam
        {params.SAMTOOLS} sort -T "{params.sort_tmp}" -O sam -n {output.fname_sam} -o {output.fname_sam_sort} >> {log.outfile} 2> >(tee -a {log.errfile} >&2)
        ## run script
        {params.PAIRED_END_READ_MERGER} {output.fname_sam_sort} {output.fname_sam_merged} {output.fname_sam_nonmerged} {input.fname_ref} >> {log.outfile} 2> >(tee -a {log.errfile} >&2)
        touch {output.fname_sam_nonmerged}
        ## sort
        rm -f '{params.sort_tmp}'.[0-9]*.bam
        {params.SAMTOOLS} sort -T "{params.sort_tmp}" -o "{output.fname_bam_merged}" "{output.fname_sam_merged}" >> {log.outfile} 2> >(tee -a {log.errfile} >&2)
        {params.SAMTOOLS} index "{output.fname_bam_merged}" >> {log.outfile} 2> >(tee -a {log.errfile} >&2)
        """


rule viloca:
    input:
        REF=(
            cohortdir("cohort_consensus.fasta")
            if config.viloca["consensus"]
            else reference_file
        ),
        BAM=(
            rules.paired_end_read_merger.output.fname_bam_merged
            if config.viloca["merge_paired_end_reads"]
            else alignment_wildcard
        ),
    output:
        SNVs="{dataset}/variants/SNVs/snvs.vcf",
        CSV="{dataset}/variants/SNVs/snv/cooccurring_mutations.csv",
        WORK_DIR=directory("{dataset}/variants/SNVs"),
    params:
        READ_LEN=read_len,
        INSERT_FILE=config.viloca["insert_bedfile"],
        MODE=config.viloca["mode"],
        SHIFT=config.viloca["shift"],
        EXTRA=config.viloca["extra"],
        VILOCA=config.applications["viloca"],
    log:
        outfile="{dataset}/variants/SNVs/viloca.out.log",
        errfile="{dataset}/variants/SNVs/viloca.err.log",
    conda:
        config.viloca["conda"]
    benchmark:
        "{dataset}/variants/SNVs/viloca.benchmark"
    threads: config.viloca["threads"]
    resources:
        disk_mb=2000,
        mem_mb=config.viloca["mem"],
        runtime=config.viloca["time"],
    shell:
        """
        let "WINDOW_SHIFTS=({params.READ_LEN} * 4/5 + {params.SHIFT}) / {params.SHIFT}"
        let "WINDOW_LEN=WINDOW_SHIFTS * {params.SHIFT}"
        echo "Windows are shifted by: ${{WINDOW_SHIFTS}} bp" > {log.outfile}
        echo "The window length is: ${{WINDOW_LEN}} bp" >> {log.outfile}

        # Get absolute path for input files
        CWD=${{PWD}}
        WORK_DIR="$(realpath -m {output.WORK_DIR})"
        BAM="$(realpath {input.BAM})"
        REF="$(realpath {input.REF})"
        OUTFILE="$(realpath -m {log.outfile})"
        ERRFILE="$(realpath -m {log.errfile})"

        # Create directory for running VILOCA
        DIR="${{WORK_DIR}}"
        if [[ ! -d "${{DIR}}" ]]; then
            echo "Creating directory ${{DIR}}" >> $OUTFILE
            mkdir -p ${{DIR}}
        fi
        # Change to the directory where VILOCA is to be executed
        cd "${{DIR}}"

        # Run VILOCA
        echo "Running VILOCA" >> $OUTFILE
        if [[ "{params.INSERT_FILE}" == "None" ]]; then
            {params.VILOCA} {params.EXTRA} -t {threads} --mode {params.MODE} -w ${{WINDOW_LEN}} -s {params.SHIFT} -b ${{BAM}} -f ${{REF}} >> $OUTFILE 2> >(tee -a $ERRFILE >&2)
        else
            INSERTFILE=${{CWD}}/{params.INSERT_FILE}
            echo "Insert file used ${{CWD}}/{params.INSERT_FILE}" >> $OUTFILE
            {params.VILOCA} {params.EXTRA} -t {threads} --mode {params.MODE} -z ${{INSERTFILE}} -b ${{BAM}} -f ${{REF}}  >> $OUTFILE 2> >(tee -a $ERRFILE >&2)
        fi

        # rename viloca output  snv/SNVs_0.010000_final.vcf --> snvs.vcf
        cp "${{WORK_DIR}}/snv/SNVs_0.010000_final.vcf" "${{WORK_DIR}}/snvs.vcf"
        """


if config.general["snv_caller"] == "shorah":

    ruleorder: snv > lofreq > viloca

elif config.general["snv_caller"] == "lofreq":

    ruleorder: lofreq > snv > viloca

elif config.general["snv_caller"] == "viloca":

    ruleorder: viloca > lofreq > snv
