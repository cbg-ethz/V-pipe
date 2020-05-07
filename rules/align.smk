import os

__author__ = "Susana Posada-Cespedes"
__author__ = "David Seifert"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"

# 1. initial consensus sequence
rule initial_vicuna:
    input:
        global_ref = reference_file,
        R1 = "{dataset}/preprocessed_data/R1.fastq",
        R2 = lambda wildcards: wildcards.dataset + "/preprocessed_data/R2.fastq" if config.input['paired'] else [],
    output:
        "{dataset}/references/vicuna_consensus.fasta"
    params:
        scratch = '1000',
        mem = config.initial_vicuna['mem'],
        time = config.initial_vicuna['time'],
        PAIRED = 'SECOND_END_FASTQ=cleaned/R2.fastq' if config.input['paired'] else '',
        PAIRED_BOOL = 'true' if config.input['paired'] else 'false',
        VICUNA = config.applications['vicuna'],
        BWA = config.applications['bwa'],
        INDELFIXER = config.applications['indelfixer'],
        CONSENSUSFIXER = config.applications['consensusfixer'],
        PICARD = config.applications['picard'],
        SAMTOOLS = config.applications['samtools'],
        WORK_DIR = "{dataset}/initial_consensus",
        FUNCTIONS = functions,
    log:
        outfile = "{dataset}/initial_consensus/vicuna.out.log",
        errfile = "{dataset}/initial_consensus/vicuna.err.log",
    conda:
        config.initial_vicuna['conda']
    benchmark:
        "{dataset}/initial_consensus/vicuna_consensus.benchmark"
    threads:
        config.initial_vicuna['threads']
    shell:
        """
        CONSENSUS_NAME={wildcards.dataset}
        CONSENSUS_NAME="${{CONSENSUS_NAME#*/}}"
        CONSENSUS_NAME="${{CONSENSUS_NAME//\//-}}"

        source {params.FUNCTIONS}

        ERRFILE=$(basename {log.errfile})
        OUTFILE=$(basename {log.outfile})

        # 1. copy initial reference for bwa
        rm -rf {params.WORK_DIR}/
        mkdir -p {params.WORK_DIR}/
        cp {input.global_ref} {params.WORK_DIR}/consensus.fasta
        cd {params.WORK_DIR}

        # 2. create bwa index
        {params.BWA} index consensus.fasta 2> >(tee $ERRFILE >&2)

        # 3. create initial alignment
        if [[ {params.PAIRED_BOOL} == "true" ]]; then
            {params.BWA} mem -t {threads} consensus.fasta ../preprocessed_data/R{{1,2}}.fastq > first_aln.sam 2> >(tee -a $ERRFILE >&2)
        else
            {params.BWA} mem -t {threads} consensus.fasta ../preprocessed_data/R1.fastq > first_aln.sam 2> >(tee -a $ERRFILE >&2)
        fi
        rm consensus.fasta.*

        # 4. remove unmapped reads
        {params.SAMTOOLS} view -b -F 4 first_aln.sam > mapped.bam 2> >(tee -a $ERRFILE >&2)
        rm first_aln.sam

        # 5. extract reads
        mkdir -p cleaned
        SamToFastq {params.PICARD} I=mapped.bam FASTQ=cleaned/R1.fastq {params.PAIRED} VALIDATION_STRINGENCY=SILENT 2> >(tee -a $ERRFILE >&2)
        rm mapped.bam

        # 6. create config file
        # NOTE: Tabs are required below
        if [[ {params.PAIRED_BOOL} == "true" ]]; then
			cat > vicuna_config.txt <<- _EOF_
				minMSize	9
				maxOverhangSize	2
				Divergence	8
				max_read_overhang	2
				max_contig_overhang	10
				pFqDir	cleaned/
				batchSize	100000
				LibSizeLowerBound	100
				LibSizeUpperBound	800
				min_output_contig_len	1000
				outputDIR	./
			_EOF_
        else
			cat > vicuna_config.txt <<- _EOF_
				minMSize	9
				maxOverhangSize	2
				Divergence	8
				max_read_overhang	2
				max_contig_overhang	10
				npFqDir	cleaned/
				batchSize	100000
				min_output_contig_len	1000
				outputDIR	./
			_EOF_
        fi

        # 7. VICUNA
        OMP_NUM_THREADS={threads} {params.VICUNA} vicuna_config.txt > $OUTFILE 2> >(tee -a $ERRFILE >&2)
        rm vicuna_config.txt
        rm -r cleaned/

        # 8. fix broken header
        sed -e 's:>dg-\([[:digit:]]\+\)\s.*:>dg-\1:g' contig.fasta > contig_clean.fasta

        # 9. InDelFixer + ConsensusFixer to polish up consensus
        for i in {{1..3}}
        do
                mv consensus.fasta old_consensus.fasta
                indelFixer {params.INDELFIXER} -i contig_clean.fasta -g old_consensus.fasta >> $OUTFILE 2> >(tee -a $ERRFILE >&2)
                sam2bam {params.SAMTOOLS} reads.sam >> $OUTFILE 2> >(tee $ERRFILE >&2)
                consensusFixer {params.CONSENSUSFIXER} -i reads.bam -r old_consensus.fasta -mcc 1 -mic 1 -d -pluralityN 0.01 >> $OUTFILE 2> >(tee $ERRFILE >&2)
        done

        sed -i -e "s/>.*/>${{CONSENSUS_NAME}}/" consensus.fasta
        echo "" >> consensus.fasta

        # 10. finally, move into place
        mkdir -p ../references
        mv {{,../references/vicuna_}}consensus.fasta
        """

rule initial_vicuna_msa:
    input:
        vicuna_refs
    output:
        "references/initial_aln_gap_removed.fasta"
    params:
        scratch = '1250',
        mem = config.initial_vicuna_msa['mem'],
        time = config.initial_vicuna_msa['time'],
        MAFFT = config.applications['mafft'],
        REMOVE_GAPS = config.applications['remove_gaps_msa'],
    log:
        outfile = "references/MAFFT_initial_aln.out.log",
        errfile = "references/MAFFT_initial_aln.err.log",
    conda:
        config.initial_vicuna_msa['conda']
    benchmark:
        "references/MAFFT_initial_aln.benchmark"
    threads:
        config.initial_vicuna_msa['threads']
    shell:
        """
        cat {input} > initial_ALL.fasta
        {params.MAFFT} --nuc --preservecase --maxiterate 1000 --localpair --thread {threads} initial_ALL.fasta > references/initial_aln.fasta 2> >(tee {log.errfile} >&2)
        rm initial_ALL.fasta

        {params.REMOVE_GAPS} references/initial_aln.fasta -o {output} -p 0.5 > {log.outfile} 2> >(tee -a {log.errfile} >&2)
        """

localrules:
    create_vicuna_initial
rule create_vicuna_initial:
    input:
        "references/initial_aln_gap_removed.fasta"
    output:
        "{dataset}/references/initial_consensus.fasta"
    params:
        EXTRACT_SEQ = config.applications['extract_seq'],
    conda:
        config.create_vicuna_initial['conda']
    shell:
        """
        CONSENSUS_NAME={wildcards.dataset}
        CONSENSUS_NAME="${{CONSENSUS_NAME#*/}}"
        CONSENSUS_NAME="${{CONSENSUS_NAME//\//-}}"

        mkdir -p {wildcards.dataset}/references/
        {params.EXTRACT_SEQ} {input} -o {output} -s "${{CONSENSUS_NAME}}"
        """

localrules:
    create_simple_initial
rule create_simple_initial:
    input:
        "references/cohort_consensus.fasta"
    output:
        "{dataset}/references/initial_consensus.fasta"
    shell:
        """
        CONSENSUS_NAME={wildcards.dataset}
        CONSENSUS_NAME="${{CONSENSUS_NAME#*/}}"
        CONSENSUS_NAME="${{CONSENSUS_NAME//\//-}}"

        mkdir -p {wildcards.dataset}/references/
        cp {input} {output}
        sed -i -e "s/>.*/>${{CONSENSUS_NAME}}/" {output}
        """

localrules:
    create_denovo_initial
rule create_denovo_initial:
    input:
        "{dataset}/references/denovo_consensus.fasta"
    output:
        "{dataset}/references/initial_consensus.fasta"
    shell:
        """
        CONSENSUS_NAME={wildcards.dataset}
        CONSENSUS_NAME="${{CONSENSUS_NAME#*/}}"
        CONSENSUS_NAME="${{CONSENSUS_NAME//\//-}}"

        mkdir -p {wildcards.dataset}/references/
        cp {input} {output}
        sed -i -e "s/>.*/>${{CONSENSUS_NAME}}/" {output}
        """


# change this to switch between VICUNA and creating a simple initial
# initial reference
ruleorder:
    create_denovo_initial > create_simple_initial > create_vicuna_initial
# ruleorder: create_vicuna_initial > create_simple_initial


# 2. aligning
def input_align(wildcards):
    list_output = []
    list_output.append(wildcards.dataset + "/preprocessed_data/R1.fastq")
    if config.input['paired']:
        list_output.append(wildcards.dataset + "/preprocessed_data/R2.fastq")
    return(list_output)


rule hmm_align:
    input:
        initial_ref = "{dataset}/references/initial_consensus.fasta",
        FASTQ = input_align,
    output:
        good_aln = temp("{dataset}/alignments/full_aln.sam"),
        reject_aln = temp("{dataset}/alignments/rejects.sam"),
        REF_ambig = "{dataset}/references/ref_ambig.fasta",
        REF_majority = "{dataset}/references/ref_majority.fasta",
    params:
        scratch = '1250',
        mem = config.hmm_align['mem'],
        time = config.hmm_align['time'],
        LEAVE_TEMP = '-l' if config.hmm_align['leave_msa_temp'] else '',
        EXTRA = config.hmm_align['extra'],
        MAFFT = config.applications['mafft'],
        NGSHMMALIGN = config.applications['ngshmmalign'],
    log:
        outfile = "{dataset}/alignments/ngshmmalign.out.log",
        errfile = "{dataset}/alignments/ngshmmalign.err.log",
    conda:
        config.hmm_align['conda']
    benchmark:
        "{dataset}/alignments/ngshmmalign.benchmark"
    threads:
        config.hmm_align['threads']
    shell:
        """
        CONSENSUS_NAME={wildcards.dataset}
        CONSENSUS_NAME="${{CONSENSUS_NAME#*/}}"
        CONSENSUS_NAME="${{CONSENSUS_NAME//\//-}}"

        # 1. clean previous run
        rm -rf   {wildcards.dataset}/alignments
        rm -f    {wildcards.dataset}/references/ref_ambig.fasta
        rm -f    {wildcards.dataset}/references/ref_majority.fasta
        mkdir -p {wildcards.dataset}/alignments
        mkdir -p {wildcards.dataset}/references

        # 2. perform alignment # -l = leave temps
        {params.NGSHMMALIGN} -v {params.EXTRA} -R {input.initial_ref} -o {output.good_aln} -w {output.reject_aln} -t {threads} -N "${{CONSENSUS_NAME}}" {params.LEAVE_TEMP} {input.FASTQ} > {log.outfile} 2> >(tee {log.errfile} >&2)

        # 3. move references into place
        mv {wildcards.dataset}/{{alignments,references}}/ref_ambig.fasta
        mv {wildcards.dataset}/{{alignments,references}}/ref_majority.fasta
        """

rule sam2bam:
    input:
        "{file}.sam"
    output:
        BAM = "{file}.bam",
        BAI = "{file}.bam.bai"
    params:
        scratch = '1250',
        mem = config.sam2bam['mem'],
        time = config.sam2bam['time'],
        SAMTOOLS = config.applications['samtools'],
        FUNCTIONS = functions,
    log:
        outfile = "{file}_sam2bam.out.log",
        errfile = "{file}_sam2bam.err.log",
    conda:
        config.sam2bam['conda']
    benchmark:
        "{file}_sam2bam.benchmark"
    threads:
        1
    shell:
        """
        # convert sam -> bam
        source {params.FUNCTIONS}
        sam2bam {params.SAMTOOLS} {input} > {log.outfile} 2> >(tee {log.errfile} >&2)
        """


# 3. construct MSA from all patient files
def construct_msa_input_files(wildcards):
    output_list = ["{}{}.fasta".format(s, wildcards.kind)
                   for s in references]
    output_list.append(reference_file)

    return output_list


rule msa:
    input:
        construct_msa_input_files
    output:
        "references/ALL_aln_{kind}.fasta"
    params:
        scratch = '1250',
        mem = config.msa['mem'],
        time = config.msa['time'],
        MAFFT = config.applications['mafft'],
    log:
        outfile = "references/MAFFT_{kind}_cohort.out.log",
        errfile = "references/MAFFT_{kind}_cohort.err.log",
    conda:
        config.msa['conda']
    benchmark:
        "references/MAFFT_{kind}_cohort.benchmark"
    threads:
        config.msa['threads']
    shell:
        """
        cat {input} > ALL_{wildcards.kind}.fasta
        {params.MAFFT} --nuc --preservecase --maxiterate 1000 --localpair --thread {threads} ALL_{wildcards.kind}.fasta > {output} 2> >(tee {log.errfile} >&2)
        rm ALL_{wildcards.kind}.fasta
        """


# 4. convert alignments to REF alignment
def get_reference_name(wildcards):
    with open(reference_file, 'r') as infile:
        reference_name = infile.readline().rstrip()
    reference_name = reference_name.split('>')[1]
    reference_name = reference_name.split(' ')[0]
    return reference_name


rule convert_to_ref:
    input:
        REF_ambig = "references/ALL_aln_ambig.fasta",
        REF_majority = "references/ALL_aln_majority.fasta",
        BAM = "{dataset}/alignments/full_aln.bam",
        REJECTS_BAM = "{dataset}/alignments/rejects.bam",
    output:
        "{dataset}/alignments/REF_aln.bam"
    params:
        scratch = '1250',
        mem = config.convert_to_ref['mem'],
        time = config.convert_to_ref['time'],
        REF_NAME = reference_name if reference_name else get_reference_name,
        CONVERT_REFERENCE = config.applications['convert_reference'],
    log:
        outfile = "{dataset}/alignments/convert_to_ref.out.log",
        errfile = "{dataset}/alignments/convert_to_ref.err.log",
    conda:
        config.convert_to_ref['conda']
    benchmark:
        "{dataset}/alignments/convert_to_ref.benchmark"
    threads:
        1
    shadow:
        "shallow"
    shell:
        """
        {params.CONVERT_REFERENCE} -t {params.REF_NAME} -m {input.REF_ambig} -i {input.BAM} -o {output} > {log.outfile} 2> >(tee {log.errfile} >&2)
        """


# 2-4. Alternative: align reads using bwa or bowtie
if config.general["aligner"] == "bwa":
    rule ref_bwa_index:
        input:
            reference_file
        output:
            "{}.bwt".format(reference_file)
        params:
            scratch = '1250',
            mem = config.ref_bwa_index['mem'],
            time = config.ref_bwa_index['time'],
            BWA = config.applications['bwa']
        log:
            outfile = "references/bwa_index.out.log",
            errfile = "references/bwa_index.err.log",
        conda:
            config.ref_bwa_index['conda']
        benchmark:
            "references/ref_bwa_index.benchmark"
        shell:
            """
            {params.BWA} index {input} 2> >(tee {log.errfile} >&2)
            """

    rule bwa_align:
        input:
            FASTQ = input_align,
            REF = reference_file,
            INDEX = "{}.bwt".format(reference_file)
        output:
            temp("{dataset}/alignments/REF_aln.sam"),
        params:
            scratch = '1250',
            mem = config.bwa_align['mem'],
            time = config.bwa_align['time'],
            EXTRA = config.bwa_align['extra'],
            FILTER = '-f 2' if config.input['paired'] else '-F 4',
            TMP_SAM = "{dataset}/alignments/tmp_aln.sam",
            BWA = config.applications['bwa'],
            SAMTOOLS = config.applications['samtools'],
        log:
            outfile = "{dataset}/alignments/bwa_align.out.log",
            errfile = "{dataset}/alignments/bwa_align.err.log",
        conda:
            config.bwa_align['conda']
        benchmark:
            "{dataset}/alignments/bwa_align.benchmark"
        threads:
            config.bwa_align['threads']
        shell:
            """
            {params.BWA} mem -t {threads} {params.EXTRA} {input.REF} {input.FASTQ} > {params.TMP_SAM} 2> >(tee {log.errfile} >&2)
            # Filter alignments: (1) remove unmapped reads (single-end) or keep only reads mapped in proper pairs (paired-end), (2) remove supplementary aligments
            {params.SAMTOOLS} view -h {params.FILTER} -F 2048 {params.TMP_SAM} > {output} 2> >(tee -a {log.errfile} >&2)
            rm {params.TMP_SAM}
            """
elif config.general["aligner"] == "bowtie":
    rule ref_bowtie_index:
        input:
            reference_file
        output:
            INDEX1 = "{}.1.bt2".format(reference_file),
            INDEX2 = "{}.2.bt2".format(reference_file),
            INDEX3 = "{}.3.bt2".format(reference_file),
            INDEX4 = "{}.4.bt2".format(reference_file),
            INDEX5 = "{}.rev.1.bt2".format(reference_file),
            INDEX6 = "{}.rev.2.bt2".format(reference_file),
        params:
            scratch = '1250',
            mem = config.ref_bowtie_index['mem'],
            time = config.ref_bowtie_index['time'],
            BOWTIE = config.applications['bowtie_idx']
        log:
            outfile = "references/bowtie_index.out.log",
            errfile = "references/bowtie_index.err.log",
        conda:
            config.ref_bowtie_index['conda']
        benchmark:
            "references/ref_bowtie_index.benchmark"
        shell:
            """
            {params.BOWTIE} {input} {input} 2> >(tee {log.errfile} >&2)
            """

    if config.input['paired']:
        rule bowtie_align:
            input:
                R1 = "{dataset}/preprocessed_data/R1.fastq.gz",
                R2 = "{dataset}/preprocessed_data/R2.fastq.gz",
                REF = reference_file,
                INDEX1 = "{}.1.bt2".format(reference_file),
                INDEX2 = "{}.2.bt2".format(reference_file),
                INDEX3 = "{}.3.bt2".format(reference_file),
                INDEX4 = "{}.4.bt2".format(reference_file),
                INDEX5 = "{}.rev.1.bt2".format(reference_file),
                INDEX6 = "{}.rev.2.bt2".format(reference_file)
            output:
                temp("{dataset}/alignments/REF_aln.sam"),
            params:
                scratch = '1250',
                mem = config.bowtie_align['mem'],
                time = config.bowtie_align['time'],
                TMP_SAM = "{dataset}/alignments/tmp_aln.sam",
                PHRED = config.bowtie_align['phred'],
                PRESET = config.bowtie_align['preset'],
                MAXINS = get_maxins, 
                EXTRA = config.bowtie_align['extra'],
                BOWTIE = config.applications['bowtie'],
                SAMTOOLS = config.applications['samtools'],
            log:
                outfile = "{dataset}/alignments/bowtie_align.out.log",
                errfile = "{dataset}/alignments/bowtie_align.err.log",
            conda:
                config.bowtie_align['conda']
            benchmark:
                "{dataset}/alignments/bowtie_align.benchmark"
            threads:
                config.bowtie_align['threads']
            shell:
                """
                {params.BOWTIE} -x {input.REF} -1 {input.R1} -2 {input.R2} {params.PHRED} {params.PRESET} -X {params.MAXINS} {params.EXTRA} -p {threads} -S {params.TMP_SAM} 2> >(tee {log.errfile} >&2)
                # Filter alignments: (1) keep only reads mapped in proper pairs, and (2) remove supplementary aligments
                {params.SAMTOOLS} view -h -f 2 -F 2048 {params.TMP_SAM} > {output} 2> >(tee -a {log.errfile} >&2)
                rm {params.TMP_SAM}
                """
    else:
        rule bowtie_align_se:
            input:
                R1 = "{dataset}/preprocessed_data/R1.fastq.gz",
                REF = reference_file,
                INDEX1 = "{}.1.bt2".format(reference_file),
                INDEX2 = "{}.2.bt2".format(reference_file),
                INDEX3 = "{}.3.bt2".format(reference_file),
                INDEX4 = "{}.4.bt2".format(reference_file),
                INDEX5 = "{}.rev.1.bt2".format(reference_file),
                INDEX6 = "{}.rev.2.bt2".format(reference_file)
            output:
                temp("{dataset}/alignments/REF_aln.sam"),
            params:
                scratch = '1250',
                mem = config.bowtie_align['mem'],
                time = config.bowtie_align['time'],
                TMP_SAM = "{dataset}/alignments/tmp_aln.sam",
                PHRED = config.bowtie_align['phred'],
                PRESET = config.bowtie_align['preset'],
                EXTRA = config.bowtie_align['extra'],
                BOWTIE = config.applications['bowtie'],
                SAMTOOLS = config.applications['samtools'],
            log:
                outfile = "{dataset}/alignments/bowtie_align.out.log",
                errfile = "{dataset}/alignments/bowtie_align.err.log",
            conda:
                config.bowtie_align['conda']
            benchmark:
                "{dataset}/alignments/bowtie_align.benchmark"
            threads:
                config.bowtie_align['threads']
            shell:
                """
                {params.BOWTIE} -x {input.REF} -U {input.R1} {params.PHRED} {params.PRESET} {params.EXTRA} -p {threads} -S {params.TMP_SAM} 2> >(tee {log.errfile} >&2)
                # Filter alignments: (1) remove unmapped reads, and (2) remove supplementary aligments
                {params.SAMTOOLS} view -h -F 4 -F 2048 {params.TMP_SAM} > {output} 2> >(tee -a {log.errfile} >&2)
                rm {params.TMP_SAM}
                """



rule consensus_sequences:
    input:
        BAM = "{dataset}/alignments/REF_aln.bam",
        REF = reference_file,
    output:
        REF_amb = "{dataset}/references/ref_ambig.fasta",
        REF_majority = "{dataset}/references/ref_majority.fasta",
    params:
        scratch = '1250',
        mem = config.consensus_sequences['mem'],
        time = config.consensus_sequences['time'],
        MIN_COVERAGE = config.consensus_sequences['min_coverage'],
        N_COVERAGE = config.consensus_sequences['n_coverage'],
        QUAL_THRD = config.consensus_sequences['qual_thrd'],
        MIN_FREQ = config.consensus_sequences['min_freq'],
        OUTDIR = "{dataset}/references",
        EXTRACT_CONSENSUS = config.applications['extract_consensus'],
    log:
        outfile = "{dataset}/references/consensus_sequences.out.log",
        errfile = "{dataset}/references/consensus_sequences.err.log",
    conda:
        config.consensus_sequences['conda']
    benchmark:
        "{dataset}/alignments/consensus.benchmark"
    threads:
        1
    shell:
        """
        CONSENSUS_NAME={wildcards.dataset}
        CONSENSUS_NAME="${{CONSENSUS_NAME#*/}}"
        CONSENSUS_NAME="${{CONSENSUS_NAME//\//-}}"

        {params.EXTRACT_CONSENSUS} -i {input.BAM} -f {input.REF} -c {params.MIN_COVERAGE} -n {params.N_COVERAGE} -q {params.QUAL_THRD} -a {params.MIN_FREQ} -N "${{CONSENSUS_NAME}}" -o {params.OUTDIR}
        """


if config.general["aligner"] == "ngshmmalign":
    ruleorder: hmm_align > consensus_sequences
    ruleorder: convert_to_ref > sam2bam
elif config.general["aligner"] == "bwa":
    ruleorder: consensus_sequences > hmm_align
    ruleorder: sam2bam > convert_to_ref
elif config.general["aligner"] == "bowtie":
    ruleorder: consensus_sequences > hmm_align
    ruleorder: sam2bam > convert_to_ref

