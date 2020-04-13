
__author__ = "Susana Posada-Cespedes"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"

# 1. Simulate haplotype sequences, including the master-haplotype sequence
rule simulate_master:
    output:
        reference_file
    params:
        scratch = '2000',
        mem = configBench.simulate_master['mem'],
        time = configBench.simulate_master['time'],
        GENOME_LEN = configBench.simulate_master['genome_length'],
        SEED = configBench.simulate_master['seed'],
        OUTDIR_HAP = "references",
        SIM_BENCH = configBench.applications['simBench'],
    log:
        outfile = "references/simulate_master.out.log",
        errfile = "references/simulate_master.out.log",
    shell:
        """
        {params.SIM_BENCH} -g {params.GENOME_LEN} -s {params.SEED} -v -oh {params.OUTDIR_HAP} -o master > >(tee {log.outfile}) 2>&1
        """

localrules:
    initial_reference
rule initial_reference:
    input:
        reference_file
    output:
        "references/cohort_consensus.fasta"
    shell:
        """
        # Get the input filename from path
        infile=$(basename {input})
        ln -s ${{infile}} {output}
        touch -h {output}
        """


def get_haplotype_seqs(wildcards):
    sample_tuple = sample_record(
        sample_name=wildcards.sample_name, date=wildcards.date)
    haplotype_seqs = sample_dict[sample_tuple]['haplotype_seqs']
    if len(haplotype_seqs) > 0 and haplotype_seqs.upper() not in ["NA", "N/A"]:
        val = haplotype_seqs
    else:
        val = ''
    return val


# def get_use_master(wildcards):
#    sample_tuple = sample_record(sample_name=wildcards.sample_name, date=wildcards.date)
#    use_master = sample_dict[sample_tuple]['use_master']
#    return '-u' if use_master else ''


def get_no_FR(wildcards):
    sample_tuple = sample_record(
        sample_name=wildcards.sample_name, date=wildcards.date)
    no_FR = sample_dict[sample_tuple]['no_frameshifts']
    return '-fr' if no_FR else ''


def get_del_len(wildcards):
    sample_tuple = sample_record(
        sample_name=wildcards.sample_name, date=wildcards.date)
    del_len = sample_dict[sample_tuple]['del_len']
    del_len = del_len.strip()
    if not del_len or del_len.upper() in ["NA", "N/A"]:
        val = ''
    else:
        val = '-dl {}'.format(del_len)
    return val


rule simulate_haplotypes:
    input:
        reference_file
    output:
        HAP = "{sample_dir}/{sample_name}/{date}/references/haplotypes/haplotypes.fasta",
    params:
        scratch = '2000',
        mem = configBench.simulate_haplotypes['mem'],
        time = configBench.simulate_haplotypes['time'],
        USE_MASTER = '-u' if configBench.simulate_haplotypes['use_master'] else '',
        HAPLOTYPE_SEQS = get_haplotype_seqs,
        NUM_HAPLOTYPES = lambda wildcards: sample_dict[sample_record(
            sample_name=wildcards.sample_name, date=wildcards.date)]['num_haplotypes'],
        MUT_RATE = lambda wildcards: sample_dict[sample_record(
            sample_name=wildcards.sample_name, date=wildcards.date)]['mut_rate'],
        DEL_RATE = lambda wildcards: sample_dict[sample_record(
            sample_name=wildcards.sample_name, date=wildcards.date)]['del_rate'],
        INS_RATE = lambda wildcards: sample_dict[sample_record(
            sample_name=wildcards.sample_name, date=wildcards.date)]['ins_rate'],
        NO_FR = get_no_FR,
        DEL_LEN = get_del_len,
        SEED = lambda wildcards: sample_dict[sample_record(
            sample_name=wildcards.sample_name, date=wildcards.date)]['seed'],
        OUTDIR_HAP = "{sample_dir}/{sample_name}/{date}/references/haplotypes",
        SIM_BENCH = configBench.applications['simBench'],
    log:
        outfile = "{sample_dir}/{sample_name}/{date}/references/haplotypes/simulate_haplotypes.out.log",
        errfile = "{sample_dir}/{sample_name}/{date}/references/haplotypes/simulate_haplotypes.out.log",
    shell:
        """
        if [[ -z {params.HAPLOTYPE_SEQS} ]]; then
            {params.SIM_BENCH} -f {input} {params.USE_MASTER} -n {params.NUM_HAPLOTYPES} -mr {params.MUT_RATE} -dr {params.DEL_RATE} -ir {params.INS_RATE} {params.NO_FR} {params.DEL_LEN} -s {params.SEED} -v -oh {params.OUTDIR_HAP} -o haplotypes > >(tee {log.outfile}) 2>&1
        else
            {params.SIM_BENCH} -f {params.HAPLOTYPE_SEQS} {params.USE_MASTER} -n {params.NUM_HAPLOTYPES} -mr {params.MUT_RATE} -dr {params.DEL_RATE} -ir {params.INS_RATE} {params.NO_FR} {params.DEL_LEN} -s {params.SEED} -v -oh {params.OUTDIR_HAP} -o haplotypes > >(tee {log.outfile}) 2>&1
        fi
        """

# 2. Simulate reads

def get_freq_params(wildcards):
    sample_tuple = sample_record(
        sample_name=wildcards.sample_name, date=wildcards.date)
    freq_dstr = sample_dict[sample_tuple]['freq_dstr']
    freq_param = sample_dict[sample_tuple]['freq_param']
    if freq_dstr == 'geom':
        val = '-gr {}'.format(freq_param)
    elif freq_dstr == 'dirichlet':
        val = '-dc {}'.format(freq_param)
    else:
        val = ''
    return val


if configBench.input['paired']:
    rule simulate_reads:
        input:
            "{sample_dir}/{sample_name}/{date}/references/haplotypes/haplotypes.fasta",
        output:
            R1_raw = "{sample_dir}/{sample_name}/{date}/raw_data/simreads_R1.fastq",
            R2_raw = "{sample_dir}/{sample_name}/{date}/raw_data/simreads_R2.fastq",
            R1 = temp(
                "{sample_dir}/{sample_name}/{date}/extracted_data/R1.fastq"),
            R2 = temp(
                "{sample_dir}/{sample_name}/{date}/extracted_data/R2.fastq"),
        params:
            scratch = '2000',
            mem = configBench.simulate_reads['mem'],
            time = configBench.simulate_reads['time'],
            NUM_HAPLOTYPES = lambda wildcards: sample_dict[sample_record(
                sample_name=wildcards.sample_name, date=wildcards.date)]['num_haplotypes'],
            COVERAGE = lambda wildcards: sample_dict[sample_record(
                sample_name=wildcards.sample_name, date=wildcards.date)]['coverage'],
            NUM_READS = '-nr' if configBench.simulate_reads['num_reads'] else '',
            HIGH_QUAL = '-q' if configBench.simulate_reads['high_quality'] else '',
            READ_LEN = lambda wildcards: sample_dict[sample_record(
                sample_name=wildcards.sample_name, date=wildcards.date)]['read_len'],
            PAIRED = '-p',
            FRAGMENT_SIZE = lambda wildcards: sample_dict[sample_record(
                sample_name=wildcards.sample_name, date=wildcards.date)]['fragment_size'],
            FREQ_DSTR = lambda wildcards: sample_dict[sample_record(
                sample_name=wildcards.sample_name, date=wildcards.date)]['freq_dstr'],
            FREQ_PARAMS = get_freq_params,
            SEED = lambda wildcards: sample_dict[sample_record(
                sample_name=wildcards.sample_name, date=wildcards.date)]['seed'],
            OUTDIR_HAP = "{sample_dir}/{sample_name}/{date}/references/haplotypes",
            OUTDIR_READS = "{sample_dir}/{sample_name}/{date}/raw_data",
            ART = configBench.applications['art'],
            SIM_BENCH = configBench.applications['simBench'],
        log:
            outfile = "{sample_dir}/{sample_name}/{date}/raw_data/simBench.out.log",
            errfile = "{sample_dir}/{sample_name}/{date}/raw_data/simBench.out.log",
        shell:
            """
            {params.SIM_BENCH} -n {params.NUM_HAPLOTYPES} -c {params.COVERAGE} {params.NUM_READS} -l {params.READ_LEN} {params.PAIRED} -m {params.FRAGMENT_SIZE} -d {params.FREQ_DSTR} {params.FREQ_PARAMS} {params.HIGH_QUAL} -art {params.ART} -s {params.SEED} -v -oh {params.OUTDIR_HAP} -or {params.OUTDIR_READS} -o reads > >(tee {log.outfile}) 2>&1

            # Move intermediate results
            mkdir -p {params.OUTDIR_READS}/reads
            mv {params.OUTDIR_READS}/*.sam {params.OUTDIR_READS}/reads
            mv {params.OUTDIR_READS}/*.aln {params.OUTDIR_READS}/reads

            # To avoid rule ambiguities: Copy simulated data sets to extracted_data folder
            cp {output.R1_raw} {output.R1}
            cp {output.R2_raw} {output.R2}
            """
else:
    rule simulate_reads:
        input:
            "{sample_dir}/{sample_name}/{date}/references/haplotypes/haplotypes.fasta",
        output:
            R1_raw = "{sample_dir}/{sample_name}/{date}/raw_data/simreads_R1.fastq",
            R1 = temp(
                "{sample_dir}/{sample_name}/{date}/extracted_data/R1.fastq"),
        params:
            scratch = '2000',
            mem = configBench.simulate_reads['mem'],
            time = configBench.simulate_reads['time'],
            NUM_HAPLOTYPES = lambda wildcards: sample_dict[sample_record(
                sample_name=wildcards.sample_name, date=wildcards.date)]['num_haplotypes'],
            COVERAGE = lambda wildcards: sample_dict[sample_record(
                sample_name=wildcards.sample_name, date=wildcards.date)]['coverage'],
            NUM_READS = '-nr' if configBench.simulate_reads['num_reads'] else '',
            READ_LEN = lambda wildcards: sample_dict[sample_record(
                sample_name=wildcards.sample_name, date=wildcards.date)]['read_len'],
            FRAGMENT_SIZE = lambda wildcards: sample_dict[sample_record(
                sample_name=wildcards.sample_name, date=wildcards.date)]['fragment_size'],
            FREQ_DSTR = lambda wildcards: sample_dict[sample_record(
                sample_name=wildcards.sample_name, date=wildcards.date)]['freq_dstr'],
            GEOM_RATIO = lambda wildcards: sample_dict[sample_record(
                sample_name=wildcards.sample_name, date=wildcards.date)]['geom_ratio'],
            SEED = lambda wildcards: sample_dict[sample_record(
                sample_name=wildcards.sample_name, date=wildcards.date)]['seed'],
            OUTDIR_HAP = "{sample_dir}/{sample_name}/{date}/references/haplotypes",
            OUTDIR_READS = "{sample_dir}/{sample_name}/{date}/raw_data",
            ART = configBench.applications['art'],
            SIM_BENCH = configBench.applications['simBench'],
        log:
            outfile = "{sample_dir}/{sample_name}/{date}/raw_data/simBench.out.log",
            errfile = "{sample_dir}/{sample_name}/{date}/raw_data/simBench.out.log",
        shell:
            """
            {params.SIM_BENCH} -n {params.NUM_HAPLOTYPES} -c {params.COVERAGE} {params.NUM_READS} -l {params.READ_LEN} -m {params.FRAGMENT_SIZE} -d {params.FREQ_DSTR} -gr {params.GEOM_RATIO} -art {params.ART} -s {params.SEED} -v -oh {params.OUTDIR_HAP} -or {params.OUTDIR_READS} -o reads > >(tee {log.outfile}) 2>&1

            # Move intermediate results
            mkdir -p {params.OUTDIR_READS}/reads
            mv {params.OUTDIR_READS}/*.sam {params.OUTDIR_READS}/reads
            mv {params.OUTDIR_READS}/*.aln {params.OUTDIR_READS}/reads

            # To avoid rule ambiguities: Copy simulated data sets to extracted_data folder
            cp {output.R1_raw} {output.R1}
            """

if configBench.general["simulate"]:
    ruleorder: simulate_reads > extract
else:
    ruleorder: extract > simulate_reads

