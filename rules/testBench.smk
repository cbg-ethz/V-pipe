
__author__ = "Susana Posada-Cespedes"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"

# 1. Test SNVs
# 1.a. Use snvs before stand-bias test


def snvfiles(wildcards):
    snv_files = []
    inferred_values = glob_wildcards(
        wildcards.dataset + "/variants/SNVs/REGION_{region}/snv/SNV.txt")
    for i in inferred_values.region:
        snv_files.append(os.path.join(wildcards.dataset, "variants", "SNVs", ''.join(
            ("REGION_", i)), "snv", "SNVs_0.010000.txt"))
    if len(snv_files) == 0:
        print("No inferred values for dataset", wildcards.dataset)
        sys.exit(1)
    return snv_files


rule aggregate_beforeSB:
    input:
        snvfiles
    output:
        TXT = temp("{dataset}/variants/SNVs/SNVs_beforeSB.txt"),
        CSV = "{dataset}/variants/SNVs/SNVs_beforeSB.csv"
    params:
        scratch = '1250',
        mem = '2000',
        time = '20'
    log:
        outfile = "{dataset}/variants/SNVs/aggregate_beforeSB.out.log",
        errfile = "{dataset}/variants/SNVs/aggregate_beforeSB.err.log"
    shell:
        """
        array=( {input} )
        num_files=( ${{#array[@]}} )
        cat {input} | sort -nk2 > {output.TXT}
        cat {output.TXT} | tr '\t' ',' > {output.CSV}
        sed -i 1i"Chromosome,Pos,Ref,Var,Frq1,Frq2,Frq3,Pst1,Pst2,Pst3,Fvar,Rvar,Ftot,Rtot,Pval" {output.CSV}
        """

# TODO wl, ws should be added
# TODO Switch between snvs and snvs_beforeSB
rule test_snv:
    input:
        # "{dataset}/variants/SNVs/SNVs_beforeSB.csv",
        SNVs = input_snv,
        HAPLOTYPE_SEQS = "{sample_dir}/{sample_name}/{date}/references/haplotypes/haplotypes.fasta",
        REF = "variants/cohort_consensus.fasta",
        REF_ALN = reference_file,
        TSV = "{sample_dir}/{sample_name}/{date}/variants/coverage_intervals.tsv",
    output:
        SNVs = "{sample_dir}/{sample_name}/{date}/variants/SNVs/true_snvs.tsv",
        PERFORMANCE = temp(
            "{sample_dir}/{sample_name}/{date}/variants/SNVs/performance.tsv")
    params:
        scratch = '2000',
        mem = config.test_snv['mem'],
        time = config.test_snv['time'],
        RE_MSA = 'true' if config.test_snv['re_msa'] else 'false',
        HAPLOTYPE_SEQS_AUX = "{sample_dir}/{sample_name}/{date}/references/haplotypes/haplotypes_aux.fasta",
        FREQ_DSTR = lambda wildcards: sample_dict[sample_record(sample_name=wildcards.sample_name, date=wildcards.date)]['freq_dstr'],
        FREQ_PARAMS = get_freq_aux,
        CALLER = '-c' if config.general['snv_caller'] != 'shorah' else '',
        OUTDIR = "{sample_dir}/{sample_name}/{date}/variants/SNVs",
        ID = lambda wildcards: f'{wildcards.sample_name}-{wildcards.date}',
        MAFFT = config.applications['mafft'],
        TEST_BENCH = config.applications['testBench'],
    log:
        outfile = "{sample_dir}/{sample_name}/{date}/variants/SNVs/testBench.out.log",
        errfile = "{sample_dir}/{sample_name}/{date}/variants/SNVs/testBench.out.log",
    conda:
        config.test_snv['conda']
    shell:
        """
        region=`tr '\n' ',' < {input.TSV}`
        if [[ -n ${{region}} ]]; then
            echo "Region(s): ${{region}}" >> {log.outfile}
            if [[ {params.RE_MSA} == "true" ]]; then
                # remove indels
                sed -e 's/-//g' {input.HAPLOTYPE_SEQS} > {params.HAPLOTYPE_SEQS_AUX}
                {params.TEST_BENCH} -f {params.HAPLOTYPE_SEQS_AUX} -s {input.SNVs} -m {input.REF} --ref {input.REF_ALN} -d {params.FREQ_DSTR} {params.FREQ_PARAMS} -r ${{region}} {params.CALLER} -t -ms -mafft {params.MAFFT} -N {params.ID} -of {output.PERFORMANCE} -od {params.OUTDIR} > >(tee -a {log.outfile}) 2>&1
            else
                {params.TEST_BENCH} -f {input.HAPLOTYPE_SEQS} -s {input.SNVs} -m {input.REF} --ref {input.REF_ALN} -d {params.FREQ_DSTR} {params.FREQ_PARAMS} -r ${{region}} {params.CALLER} -t -N {params.ID} -of {output.PERFORMANCE} -od {params.OUTDIR} > >(tee -a {log.outfile}) 2>&1
            fi
        else
            echo "No called SNVs"
            touch {output}
        fi
        """

rule aggregate:
    input:
        expand("{dataset}/variants/SNVs/performance.tsv", dataset=datasets)
    output:
        "variants/SNV_calling_performance.tsv"
    params:
        scratch = '1250',
        mem = config.aggregate['mem'],
        time = config.aggregate['time']
    log:
        outfile = "variants/SNV_calling_performance.out.log",
        errfile = "variants/SNV_calling_performance.err.log"
    shell:
        """
        awk FNR!=1 {input} > {output}
        sed -i 1i"ID\tTP\tFP\tFN\tTN" {output}
        """
