__author__ = "Susana Posada-Cespedes"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"

# 1. Test SNVs
# 1.a. Use snvs before stand-bias test


def snvfiles(wildcards):
    snv_files = []
    inferred_values = glob_wildcards(
        wildcards.dataset + "/variants/SNVs/REGION_{region}/snv/SNV.txt"
    )
    for i in inferred_values.region:
        snv_files.append(
            os.path.join(
                wildcards.dataset,
                "variants",
                "SNVs",
                "".join(("REGION_", i)),
                "snv",
                "SNVs_0.010000.txt",
            )
        )
    if len(snv_files) == 0:
        print("No inferred values for dataset", wildcards.dataset)
        sys.exit(1)
    return snv_files


rule alignment_bias:
    input:
        REF=reference_file,
        BAM="{sample_dir}/{sample_name}/{date}/alignments/REF_aln.bam",
        R1gz="{sample_dir}/{sample_name}/{date}/preprocessed_data/R1.fastq.gz",
        HAPLOTYPE_SEQS=(
            "{sample_dir}/{sample_name}/{date}/references/haplotypes/haplotypes.fasta"
        ),
    output:
        "{sample_dir}/{sample_name}/{date}/alignments/alignment_bias.tsv",
    params:
        scratch="2000",
        mem=config.alignment_bias["mem"],
        time=config.alignment_bias["time"],
        PAIRED="-p" if config.input["paired"] else "",
        ID=lambda wildcards: f"{wildcards.sample_name}-{wildcards.date}",
        ALIGNMENT_BIAS=config.applications["alignmentBias"],
    log:
        outfile="{sample_dir}/{sample_name}/{date}/alignments/alignment_bias.out.log",
        errfile="{sample_dir}/{sample_name}/{date}/alignments/alignment_bias.out.log",
    conda:
        config.alignment_bias["conda"]
    threads: 1
    shell:
        """
        {params.ALIGNMENT_BIAS} -r {input.REF} -b {input.BAM} -f <(zcat {input.R1gz}) --hap {input.HAPLOTYPE_SEQS} {params.PAIRED} -N {params.ID} -o {output}
        """


rule aggregate_alignment_bias:
    input:
        expand("{dataset}/alignments/alignment_bias.tsv", dataset=datasets),
    output:
        "stats/alignment_bias.tsv",
    params:
        scratch="1250",
        mem=config.aggregate["mem"],
        time=config.aggregate["time"],
    log:
        outfile="stats/alignment_bias.out.log",
        errfile="stats/alignment_bias.out.log",
    shell:
        """
        awk FNR!=1 {input} > {output}
        sed -i 1i"SampleID\tHaplotypeID\tDivergence\tPercent-aligned\tPercent-bases-aligne\n" {output}
        """


rule aggregate_beforeSB:
    input:
        snvfiles,
    output:
        TXT=temp("{dataset}/variants/SNVs/SNVs_beforeSB.txt"),
        CSV="{dataset}/variants/SNVs/SNVs_beforeSB.csv",
    params:
        scratch="1250",
        mem="2000",
        time="20",
    log:
        outfile="{dataset}/variants/SNVs/aggregate_beforeSB.out.log",
        errfile="{dataset}/variants/SNVs/aggregate_beforeSB.err.log",
    shell:
        """
        array=( {input} )
        num_files=( ${{#array[@]}} )
        cat {input} | sort -nk2 > {output.TXT}
        cat {output.TXT} | tr '\t' ',' > {output.CSV}
        sed -i 1i"Chromosome,Pos,Ref,Var,Frq1,Frq2,Frq3,Pst1,Pst2,Pst3,Fvar,Rvar,Ftot,Rtot,Pval" {output.CSV}
        """


# TODO Switch between snvs and snvs_beforeSB
rule test_snv:
    input:
        INPUT=input_snv,
        HAPLOTYPE_SEQS=(
            "{sample_dir}/{sample_name}/{date}/references/haplotypes/haplotypes.fasta"
        ),
        REF=(
            "variants/cohort_consensus.fasta"
            if config.snv["consensus"]
            else reference_file
        ),
        REF_ALN=reference_file,
    output:
        temp("{sample_dir}/{sample_name}/{date}/variants/SNVs/performance.tsv"),
    params:
        scratch="2000",
        mem=config.test_snv["mem"],
        time=config.test_snv["time"],
        RE_MSA="true" if config.test_snv["re_msa"] else "false",
        HAPLOTYPE_SEQS_AUX="{sample_dir}/{sample_name}/{date}/references/haplotypes/haplotypes_aux.fasta",
        FREQ_DSTR=lambda wildcards: sample_dict[
            sample_record(sample_name=wildcards.sample_name, date=wildcards.date)
        ]["freq_dstr"],
        FREQ_PARAMS=get_freq_aux,
        CALLER_OPT=(
            "-ci" if config.general["snv_caller"] == "shorah" else "--no-shorah -cf"
        ),
        WINDOW_LEN=window_len,
        WINDOW_SHIFT=config.snv["shift"],
        OUTDIR="{sample_dir}/{sample_name}/{date}/variants/SNVs",
        ID=lambda wildcards: f"{wildcards.sample_name}-{wildcards.date}",
        MAFFT=config.applications["mafft"],
        EXTRA=config.test_snv["extra"],
        TEST_BENCH=config.applications["testBench"],
    log:
        outfile="{sample_dir}/{sample_name}/{date}/variants/SNVs/testBench.out.log",
        errfile="{sample_dir}/{sample_name}/{date}/variants/SNVs/testBench.out.log",
    conda:
        config.test_snv["conda"]
    threads: 1
    shell:
        """
        if [[ {params.RE_MSA} == "true" ]]; then
            # remove indels
            sed -e 's/-//g' {input.HAPLOTYPE_SEQS} > {params.HAPLOTYPE_SEQS_AUX}
            {params.TEST_BENCH} -f {params.HAPLOTYPE_SEQS_AUX} \
                -s {input.INPUT[0]} \
                -m {input.REF} \
                --ref {input.REF_ALN} \
                -d {params.FREQ_DSTR} \
                {params.FREQ_PARAMS} \
                {params.CALLER_OPT} {input.INPUT[1]} \
                -wl {params.WINDOW_LEN} -ws {params.WINDOW_SHIFT} \
                -t -ms -mafft {params.MAFFT} \
                -N {params.ID} \
                -of {output} \
                -od {params.OUTDIR} \
                {params.EXTRA} > >(tee -a {log.outfile}) 2>&1

            rm {params.HAPLOTYPE_SEQS_AUX}
        else
            {params.TEST_BENCH} -f {input.HAPLOTYPE_SEQS} \
                -s {input.INPUT[0]} \
                -m {input.REF} \
                --ref {input.REF_ALN} \
                -d {params.FREQ_DSTR} \
                {params.FREQ_PARAMS} \
                {params.CALLER_OPT} {input.INPUT[1]} \
                -wl {params.WINDOW_LEN} -ws {params.WINDOW_SHIFT} \
                -t \
                -N {params.ID} \
                -of {output} \
                -od {params.OUTDIR} \
                {params.EXTRA} > >(tee -a {log.outfile}) 2>&1
        fi
        """


rule compare_snv:
    input:
        REF=(
            "variants/cohort_consensus.fasta"
            if config.snv["consensus"]
            else reference_file
        ),
        REF_ALN=reference_file,
        TSV=input_tsv,
    output:
        temp("{sample_dir}/{sample_name}/{date}/variants/SNVs/{kind}/performance.tsv"),
    params:
        scratch="2000",
        mem=config.test_snv["mem"],
        time=config.test_snv["time"],
        RE_MSA="true" if config.test_snv["re_msa"] else "false",
        SNVs="{sample_dir}/{sample_name}/{date}/variants/SNVs/snvs.vcf",
        HAPLOTYPE_SEQS=(
            "{sample_dir}/{sample_name}/{date}/references/haplotypes/haplotypes.fasta"
        ),
        HAPLOTYPE_SEQS_AUX="{sample_dir}/{sample_name}/{date}/references/haplotypes/haplotypes_{kind}.fasta",
        FREQ_DSTR=lambda wildcards: sample_dict[
            sample_record(sample_name=wildcards.sample_name, date=wildcards.date)
        ]["freq_dstr"],
        FREQ_PARAMS=get_freq_aux,
        CALLER_OPT=(
            "--caller lofreq"
            if config.general["snv_caller"] == "lofreq"
            else "--no-expansion"
        ),
        OUTDIR="{sample_dir}/{sample_name}/{date}/variants/SNVs/{kind}",
        ID=lambda wildcards: f"{wildcards.sample_name}-{wildcards.date}",
        MAFFT=config.applications["mafft"],
        TEST_BENCH=config.applications["testBench"],
    log:
        outfile=(
            "{sample_dir}/{sample_name}/{date}/variants/SNVs/{kind}/testBench.out.log"
        ),
        errfile=(
            "{sample_dir}/{sample_name}/{date}/variants/SNVs/{kind}/testBench.out.log"
        ),
    conda:
        config.test_snv["conda"]
    threads: 1
    shell:
        """
        if [[ {params.RE_MSA} == "true" ]]; then
            # remove indels
            cp {params.HAPLOTYPE_SEQS} {params.OUTDIR}/hap_tmp.fasta
            sed -e 's/-//g' {params.OUTDIR}/hap_tmp.fasta > {params.HAPLOTYPE_SEQS_AUX}
            {params.TEST_BENCH} -f {params.HAPLOTYPE_SEQS_AUX} \
                -s {params.SNVs} \
                -m {input.REF} \
                --ref {input.REF_ALN} \
                -d {params.FREQ_DSTR} \
                {params.FREQ_PARAMS} \
                -ci {input.TSV} \
                {params.CALLER_OPT} -t \
                -ms -mafft {params.MAFFT} \
                -N {params.ID} \
                -of {output} \
                -od {params.OUTDIR} \
                {params.EXTRA} > >(tee -a {log.outfile}) 2>&1

            rm {params.OUTDIR}/hap_tmp.fasta
            rm {params.HAPLOTYPE_SEQS_AUX}
        else
            {params.TEST_BENCH} -f {params.HAPLOTYPE_SEQS} \
                -s {params.SNVs} \
                -m {input.REF} \
                --ref {input.REF_ALN} \
                -d {params.FREQ_DSTR} \
                {params.FREQ_PARAMS} \
                -ci {input.TSV} \
                {params.CALLER_OPT} -t \
                -N {params.ID} \
                -of {output} \
                -od {params.OUTDIR} \
                {params.EXTRA} > >(tee -a {log.outfile}) 2>&1
        fi
        """


rule aggregate:
    input:
        expand("{dataset}/variants/SNVs/performance.tsv", dataset=datasets),
    output:
        "variants/SNV_calling_performance.tsv",
    params:
        scratch="1250",
        mem=config.aggregate["mem"],
        time=config.aggregate["time"],
    log:
        outfile="variants/SNV_calling_performance.out.log",
        errfile="variants/SNV_calling_performance.out.log",
    shell:
        """
        awk FNR!=1 {input} > {output}
        sed -i 1i"ID\tTP\tFP\tFN\tTN" {output}
        """


rule aggregate_kind:
    input:
        expand("{dataset}/variants/SNVs/{{kind}}/performance.tsv", dataset=datasets),
    output:
        "variants/SNV_calling_performance_{kind}.tsv",
    params:
        scratch="1250",
        mem=config.aggregate["mem"],
        time=config.aggregate["time"],
    log:
        outfile="variants/SNV_calling_performance_{kind}.out.log",
        errfile="variants/SNV_calling_performance_{kind}.out.log",
    shell:
        """
        awk FNR!=1 {input} > {output}
        sed -i 1i"ID\tTP\tFP\tFN\tTN" {output}
        """
