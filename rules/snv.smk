import os


# 1. Call single nucleotide variants

def window_lengths(wildcards):
    window_len = []
    for p in patient_list:
        read_len = patient_dict[p]
        aux = int(
            (read_len * 4 / 5 + config.snv['shift']) / config.snv['shift'])
        window_len.append(str(aux * config.snv['shift']))

    window_len = ','.join(window_len)
    return window_len


def shifts(wildcards):
    shifts = []
    for p in patient_list:
        read_len = patient_dict[p]
        aux = int(
            (read_len * 4 / 5 + config.snv['shift']) / config.snv['shift'])
        shifts.append(str(aux))

    shifts = ','.join(shifts)
    return shifts


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
        REF = "variants/cohort_consensus.fasta",
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
        SHIFT = config.snv['shift'],
        KEEP_FILES = 'true' if config.snv['keep_files'] else 'false',
        WORK_DIR = "{dataset}/variants/SNVs",
        SHORAH = config.applications['shorah'],
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

        # Run ShoRAH in each of the predetermined regions (regions with sufficient coverage)
        LINE_COUNTER=0
        FILES=""
        FILES_VCF=""
        while read -r region || [[ -n ${{region}} ]]
        do
            echo "Running ShoRAH on region: ${{region}}" >> $OUTFILE
            LINE_COUNTER=$(( $LINE_COUNTER + 1))
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
            cd ${{DIR}}

            # NOTE: Execution command for ShoRAH2 valid from v1.99.0 and above
            {params.SHORAH} -w ${{WINDOW_LEN}} -x 100000 -r ${{region}} -R 42 -b ${{BAM}} -f ${{REF}} >> $OUTFILE 2> >(tee -a $ERRFILE >&2)
            {params.BCFTOOLS} view ${{DIR}}/snv/SNVs_0.010000_final.vcf -Oz -o ${{DIR}}/snv/SNVs_0.010000_final.vcf.gz
            {params.BCFTOOLS} index ${{DIR}}/snv/SNVs_0.010000_final.vcf.gz
            if [[ -f ${{DIR}}/snv/SNVs_0.010000_final.csv ]]; then
                FILES="$FILES ${{DIR}}/snv/SNVs_0.010000_final.csv"
                FILES_VCF="$FILES_VCF ${{DIR}}/snv/SNVs_0.010000_final.vcf.gz"
            else
                echo "ERROR: unsuccesful execution of ShoRAH" 2> >(tee -a $ERRFILE >&2)
                exit 1
            fi

            # Change back to working directory
            cd ${{CWD}}
        done < {input.TSV}

        # Aggregate results from different regions
        if [[ -z ${{FILES}} ]]; then
            if [[ ${{LINE_COUNTER}} > 0 ]]; then
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

rule snvclean:
    params:
        DIR = config.input['datadir']
    shell:
        """
        rm -rf {params.DIR}/*/*/variants/SNVs
        """

rule lofreq:
    input:
        REF = "variants/cohort_consensus.fasta",
        BAM = "{dataset}/alignments/REF_aln.bam",
        TSV = "{dataset}/variants/coverage_intervals.tsv",
    output:
        BAM = "{dataset}/variants/SNVs/REF_aln_indelqual.bam",
        SNVs = "{dataset}/variants/SNVs/snvs.vcf"
    params:
        scratch = '2000',
        mem = config.lofreq['mem'],
        time = config.lofreq['time'],
        OUTDIR = "{dataset}/variants/SNVs",
        SAMTOOLS = config.applications['samtools'],
        BCFTOOLS = config.applications['bcftools'],
        LOFREQ = config.applications['lofreq'],
    log:
        outfile = "{dataset}/variants/SNVs/lofreq.out.log",
        errfile = "{dataset}/variants/SNVs/lofreq.out.log",
    conda:
        config.lofreq['conda']
    shell:
        """
        # Add qualities to indels
        {params.LOFREQ} indelqual --dindel -f {input.REF} -o {output.BAM} --verbose {input.BAM} > >(tee -a {log.outfile}) 2>&1
        # Index bam file
        {params.SAMTOOLS} index {output.BAM}

        # Run Lofreq
        # NOTE: lofreq reads the region as a closed interval and uses 1-based indexing
        LINE_COUNTER=0
        FILES=""
        while read -r region || [[ -n ${{region}} ]]
        do
            LINE_COUNTER=$(( $LINE_COUNTER + 1 ))
            echo "Running Lofreq in region: ${{region}}" >> {log.outfile}
            OUTFILE_REGION={params.OUTDIR}/snvs_${{LINE_COUNTER}}.vcf
            {params.LOFREQ} call --call-indels -f {input.REF} -r ${{region}} -o ${{OUTFILE_REGION}} --verbose {output.BAM} > >(tee -a {log.outfile}) 2>&1
            {params.BCFTOOLS} view ${{OUTFILE_REGION}} -Oz -o ${{OUTFILE_REGION}}.gz
            {params.BCFTOOLS} index ${{OUTFILE_REGION}}.gz
            FILES="$FILES ${{OUTFILE_REGION}}.gz"
        done < {input.TSV}

        # Aggregate results from different regions
        echo "Intermediate vcf files: ${{FILES}}" >> {log.outfile}
        {params.BCFTOOLS} concat -o {params.OUTDIR}/snvs_tmp.vcf ${{FILES}}
        {params.BCFTOOLS} sort {params.OUTDIR}/snvs_tmp.vcf -o {output.SNVs}
        rm -f {params.OUTDIR}/snvs_tmp.vcf
        """

if config.general["snv_caller"] == "shorah":
    ruleorder: snv > lofreq
elif config.general["snv_caller"] == "lofreq":
    ruleorder: lofreq > snv

