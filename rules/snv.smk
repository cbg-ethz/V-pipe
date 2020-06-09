import os

__author__ = "Susana Posada-Cespedes"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"


# 1. Call single nucleotide variants
rule coverage_intervals:
    input:
        BAM = expand("{dataset}/alignments/REF_aln.bam", dataset=datasets),
        TSV = "variants/coverage.tsv",
    output:
        "variants/coverage_intervals.tsv"
    params:
        scratch = '1250',
        mem = config.coverage_intervals['mem'],
        time = config.coverage_intervals['time'],
        WINDOW_LEN = window_lengths,
        COVERAGE = config.coverage_intervals['coverage'],
        OVERLAP = '' if config.coverage_intervals['overlap'] else '-cf variants/coverage.tsv',
        SHIFT = shifts,
        NAMES = IDs,
        LIBERAL = '-e' if config.coverage_intervals['liberal'] else '',
        EXTRACT_COVERAGE_INTERVALS = config.applications['extract_coverage_intervals']
    log:
        outfile = "variants/coverage_intervals.out.log",
        errfile = "variants/coverage_intervals.out.log",
    conda:
        config.coverage_intervals['conda']
    benchmark:
        "variants/coverage_intervals.benchmark"
    threads:
        config.coverage_intervals['threads']
    shell:
        """
        {params.EXTRACT_COVERAGE_INTERVALS} -c {params.COVERAGE} -w {params.WINDOW_LEN} -s {params.SHIFT} -N {params.NAMES} {params.LIBERAL} {params.OVERLAP} -t {threads} -o {output} {input.BAM} > >(tee {log.outfile}) 2>&1
        """

localrules:
    shorah_regions
rule shorah_regions:
    input:
        "variants/coverage_intervals.tsv"
    output:
        temp(
            expand("{dataset}/variants/coverage_intervals.tsv", dataset=datasets))
    params:
        scratch = '1250',
    threads:
        1
    run:
        with open(input[0], 'r') as infile:
            for line in infile:
                parts = line.rstrip().split('\t')
                patientID = parts[0].split('-')
                sample_date = patientID[-1]
                patientID = '-'.join(patientID[:-1])
                if len(parts) == 2:
                    regions = parts[1].split(',')
                else:
                    regions = []

                with open(os.path.join(config.input['datadir'], patientID, sample_date, "variants", "coverage_intervals.tsv"), 'w') as outfile:
                    outfile.write('\n'.join(regions))


def read_len(wildcards):
    parts = wildcards.dataset.split('/')
    patient_ID = parts[1]
    date = parts[2]
    patient_tuple = patient_record(patient_id=patient_ID, date=date)
    read_len = patient_dict[patient_tuple]
    return read_len

rule snv:
    input:
        REF = "variants/cohort_consensus.fasta" if config.snv['consensus'] else reference_file,
        BAM = "{dataset}/alignments/REF_aln.bam",
        TSV = "{dataset}/variants/coverage_intervals.tsv",
    output:
        CSV = "{dataset}/variants/SNVs/snvs.csv",
        VCF = "{dataset}/variants/SNVs/snvs.vcf"
    params:
        scratch = '1250',
        mem = config.snv['mem'],
        time = config.snv['time'],
        READ_LEN = read_len,
        ALPHA = config.snv['alpha'],
        IGNORE_INDELS = '--ignore_indels' if config.snv['ignore_indels'] else '',
        COVERAGE = config.snv['coverage'],
        SHIFT = config.snv['shift'],
        KEEP_FILES = 'true' if config.snv['keep_files'] else 'false',
        WORK_DIR = "{dataset}/variants/SNVs",
        LOCALSCRATCH = config.snv['localscratch'],
        SHORAH = config.applications['shorah'],
        COVINT = config.coverage_intervals['coverage'],
        BCFTOOLS = config.applications['bcftools']
    log:
        outfile = "{dataset}/variants/SNVs/shorah.out.log",
        errfile = "{dataset}/variants/SNVs/shorah.err.log",
    conda:
        config.snv['conda']
    benchmark:
        "{dataset}/variants/SNVs/shorah.benchmark"
    threads:
        config.snv['threads']
    shell:
        """
        let "WINDOW_SHIFTS=({params.READ_LEN} * 4/5 + {params.SHIFT}) / {params.SHIFT}"
        let "WINDOW_LEN=WINDOW_SHIFTS * {params.SHIFT}"

        echo "Windows are shifted by: ${{WINDOW_SHIFTS}} bp" > {log.outfile}
        echo "The window length is: ${{WINDOW_LEN}} bp" >> {log.outfile}

        # Get absolute path for input files
        CWD=${{PWD}}
        BAM=${{PWD}}/{input.BAM}
        REF=${{PWD}}/{input.REF}
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
        FILES=""
        FILES_VCF=""
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
            {params.SHORAH} -a {params.ALPHA} -w ${{WINDOW_LEN}} -x 100000 {params.IGNORE_INDELS} -c {params.COVERAGE} -r ${{region}} -R 42 -b ${{BAM}} -f ${{REF}} >> $OUTFILE 2> >(tee -a $ERRFILE >&2)
            if [[ -n "{params.LOCALSCRATCH}" ]]; then
                # copyback from localscratch
                rsync -auq "{params.LOCALSCRATCH}/REGION_${{LINE_COUNTER}}" "${{WORK_DIR}}"
                cd ${{DIR}}
            fi
            if [[ -s ${{DIR}}/reads.fas && -f ${{DIR}}/snv/SNVs_0.010000_final.csv ]]; then
                # Non empty reads: run had enough data and should have produced SNVs
                {params.BCFTOOLS} view ${{DIR}}/snv/SNVs_0.010000_final.vcf -Oz -o ${{DIR}}/snv/SNVs_0.010000_final.vcf.gz
                {params.BCFTOOLS} index ${{DIR}}/snv/SNVs_0.010000_final.vcf.gz
                FILES="$FILES ${{DIR}}/snv/SNVs_0.010000_final.csv"
                FILES_VCF="$FILES_VCF ${{DIR}}/snv/SNVs_0.010000_final.vcf.gz"
            else
                # if we have disabled coverage intervales entirely, the first and only line might have no reads
                # (e.g.: in negative controls )
                if (( {params.COVINT} == 0 && LINE_COUNTER == 1 )) && [[ -f ${{DIR}}/reads.fas && ( ! -s ${{DIR}}/reads.fas ) ]]; then
                     echo "No reads while coverage intervals disabled (possible negative control sample)" 2> >(tee -a $ERRFILE >&2)
                     cd ${{CWD}}
                     (( --LINE_COUNTER )) || true # Strict mode : (( 0 )) = fail
                     break
                fi
                echo "ERROR: unsuccesful execution of ShoRAH" 2> >(tee -a $ERRFILE >&2)
                exit 1
            fi

            # Change back to working directory
            cd ${{CWD}}
        done < {input.TSV}

        # Aggregate results from different regions
        if [[ -z ${{FILES}} ]]; then
            if (( LINE_COUNTER > 0 )); then
                echo "ERROR: unsuccesful execution of ShoRAH" 2> >(tee -a {log.errfile} >&2)
                exit 1
            else
                echo "No alignment region reports sufficient coverage" >> {log.outfile}
                touch {output.CSV}
                touch {output.VCF}
            fi
        else
            echo "Intermediate csv files: ${{FILES}}" >> {log.outfile}
            echo "Intermediate vcf files: ${{FILES_VCF}}" >> {log.outfile}
            cat ${{FILES}} | sort -t, -nk2 | tail -n +${{LINE_COUNTER}} > {output.CSV}
            {params.BCFTOOLS} concat -o ${{WORK_DIR}}/snvs_tmp.vcf ${{FILES_VCF}}
            {params.BCFTOOLS} sort ${{WORK_DIR}}/snvs_tmp.vcf  -o {output.VCF}
            rm -f ${{WORK_DIR}}/snvs_tmp.vcf
        fi
        """

rule samtools_index:
    input:
        "{file}.fasta",
    output:
        "{file}.fasta.fai",
    params:
        scratch = '2000',
        mem = config.samtools_index['mem'],
        time = config.samtools_index['time'],
        SAMTOOLS = config.applications['samtools'],
    log:
        outfile = "{file}_samtools_index.out.log",
        errfile = "{file}_samtools_index.err.log",
    conda:
        config.samtools_index['conda']
    benchmark:
        "{file}_samtools_index.benchmark"
    shell:
        """
        {params.SAMTOOLS} faidx {input} -o {output} > {log.outfile} 2> >(tee -a {log.errfile} >&2)
        """

rule lofreq:
    input:
        REF = "variants/cohort_consensus.fasta" if config.lofreq['consensus'] else reference_file,
        REF_IDX = "variants/cohort_consensus.fasta.fai",
        BAM = "{dataset}/alignments/REF_aln.bam",
    output:
        BAM = "{dataset}/variants/SNVs/REF_aln_indelqual.bam",
        SNVs = "{dataset}/variants/SNVs/snvs.vcf"
    params:
        scratch = '2000',
        mem = config.lofreq['mem'],
        time = config.lofreq['time'],
        OUTDIR = "{dataset}/variants/SNVs",
        EXTRA = config.lofreq['extra'],
        SAMTOOLS = config.applications['samtools'],
        LOFREQ = config.applications['lofreq'],
    log:
        outfile = "{dataset}/variants/SNVs/lofreq.out.log",
        errfile = "{dataset}/variants/SNVs/lofreq.err.log",
    conda:
        config.lofreq['conda']
    benchmark:
        "{dataset}/variants/SNVs/lofreq.benchmark"
    shell:
        """
        # Add qualities to indels
        {params.LOFREQ} indelqual --dindel -f {input.REF} -o {output.BAM} --verbose {input.BAM} > {log.outfile} 2> >(tee -a {log.errfile} >&2)
        # Index bam file
        {params.SAMTOOLS} index {output.BAM} 2> >(tee {log.errfile} >&2)

        # Run Lofreq
        echo "Running LoFreq" >> {log.outfile}
        {params.LOFREQ} call {params.EXTRA} --call-indels -f {input.REF} -o {output.SNVs} --verbose {output.BAM} >> {log.outfile} 2> >(tee -a {log.errfile} >&2)
        """

if config.general["snv_caller"] == "shorah":
    ruleorder: snv > lofreq
elif config.general["snv_caller"] == "lofreq":
    ruleorder: lofreq > snv

