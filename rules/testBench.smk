
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


def get_freq_aux(wildcards):
    sample_tuple = sample_record(
        sample_name=wildcards.sample_name, date=wildcards.date)
    freq_dstr = sample_dict[sample_tuple]['freq_dstr']
    freq_param = sample_dict[sample_tuple]['freq_param']
    if freq_dstr == 'geom':
        val = '-gr {}'.format(freq_param)
    elif freq_dstr == 'dirichlet':
        infile = os.path.join(wildcards.sample_dir, wildcards.sample_name,
                              wildcards.date, "references/haplotypes/haplotype_frequencies.fasta")
        val = '-df {}'.format(infile)
    else:
        val = ''
    return val


def input_snv(wildcards):
    if config.general['snv_caller'] == 'shorah':
        output = os.path.join(wildcards.sample_dir, wildcards.sample_name,
                              wildcards.date, "variants/SNVs/snvs.csv")
    elif config.general['snv_caller'] == 'lofreq':
        output = os.path.join(wildcards.sample_dir, wildcards.sample_name,
                              wildcards.date, "variants/SNVs/snvs.vcf")
    return output

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
        "{sample_dir}/{sample_name}/{date}/variants/SNVs/true_snvs.tsv"
    params:
        scratch = '2000',
        mem = configBench.test_snv['mem'],
        time = configBench.test_snv['time'],
        RE_MSA = 'true' if configBench.test_snv['re_msa'] else 'false',
        HAPLOTYPE_SEQS_AUX = "{sample_dir}/{sample_name}/{date}/references/haplotypes/haplotypes_aux.fasta",
        FREQ_DSTR = lambda wildcards: sample_dict[sample_record(sample_name=wildcards.sample_name, date=wildcards.date)]['freq_dstr'],
        FREQ_PARAMS = get_freq_aux,
        CALLER = '-c' if configBench.general['snv_caller'] != 'shorah' else '',
        OUTDIR = "{sample_dir}/{sample_name}/{date}/variants/SNVs",
        MAFFT = configBench.applications['mafft'],
        TEST_BENCH = configBench.applications['testBench'],
    log:
        outfile = "{sample_dir}/{sample_name}/{date}/variants/SNVs/testBench.out.log",
        errfile = "{sample_dir}/{sample_name}/{date}/variants/SNVs/testBench.out.log",
    conda:
        configBench.test_snv['conda']
    shell:
        """
        region=`tr '\n' ',' < {input.TSV}`
        if [[ -n ${{region}} ]]; then
            echo "Region(s): ${{region}}" >> {log.outfile}
            if [[ {params.RE_MSA} == "true" ]]; then
                # remove indels
                sed -e 's/-//g' {input.HAPLOTYPE_SEQS} > {params.HAPLOTYPE_SEQS_AUX}
                {params.TEST_BENCH} -f {params.HAPLOTYPE_SEQS_AUX} -s {input.SNVs} -m {input.REF} --ref {input.REF_ALN} -d {params.FREQ_DSTR} {params.FREQ_PARAMS} -r ${{region}} {params.CALLER} -t -ms -mafft {params.MAFFT} -of {output} -od {params.OUTDIR} > >(tee -a {log.outfile}) 2>&1
            else
                {params.TEST_BENCH} -f {input.HAPLOTYPE_SEQS} -s {input.SNVs} -m {input.REF} --ref {input.REF_ALN} -d {params.FREQ_DSTR} {params.FREQ_PARAMS} -r ${{region}} {params.CALLER} -t -of {output} -od {params.OUTDIR} > >(tee -a {log.outfile}) 2>&1
            fi
        else
            echo "No called SNVs"
            touch {output}
        fi
        """
