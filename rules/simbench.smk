__author__ = "Susana Posada-Cespedes"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"


# 1. Simulate haplotype sequences, including the master-haplotype sequence
rule simulate_master:
    output:
        reference_file,
    params:
        scratch="2000",
        mem=config.simulate_master["mem"],
        time=config.simulate_master["time"],
        GENOME_LEN=config.simulate_master["genome_length"],
        SEED=config.simulate_master["seed"],
        OUTDIR_HAP="references",
        SIM_BENCH=config.applications["simBench"],
    log:
        outfile="references/simulate_master.out.log",
        errfile="references/simulate_master.out.log",
    conda:
        config.simulate_master["conda"]
    shell:
        """
        {params.SIM_BENCH} -g {params.GENOME_LEN} -s {params.SEED} -v -oh {params.OUTDIR_HAP} -o master > >(tee {log.outfile}) 2>&1
        """


localrules:
    initial_reference,


rule initial_reference:
    input:
        reference_file,
    output:
        "references/cohort_consensus.fasta",
    shell:
        """
        # Get the input filename from path
        infile=$(basename {input})
        ln -s ${{infile}} {output}
        touch -h {output}
        """


rule simulate_haplotypes:
    input:
        reference_file,
    output:
        HAP="{sample_dir}/{sample_name}/{date}/references/haplotypes/haplotypes.fasta",
    params:
        scratch="2000",
        mem=config.simulate_haplotypes["mem"],
        time=config.simulate_haplotypes["time"],
        HAPLOTYPE_SEQS=get_haplotype_seqs,
        NUM_HAPLOTYPES=lambda wildcards: sample_dict[
            sample_record(sample_name=wildcards.sample_name, date=wildcards.date)
        ]["num_haplotypes"],
        MUT_RATE=lambda wildcards: sample_dict[
            sample_record(sample_name=wildcards.sample_name, date=wildcards.date)
        ]["mut_rate"],
        DEL_RATE=lambda wildcards: sample_dict[
            sample_record(sample_name=wildcards.sample_name, date=wildcards.date)
        ]["del_rate"],
        INS_RATE=lambda wildcards: sample_dict[
            sample_record(sample_name=wildcards.sample_name, date=wildcards.date)
        ]["ins_rate"],
        NO_FR=get_no_FR,
        DEL_LEN=get_del_len,
        SEED=lambda wildcards: sample_dict[
            sample_record(sample_name=wildcards.sample_name, date=wildcards.date)
        ]["seed"],
        TREE_LIKE="--tree-like" if config.simulate_haplotypes["tree_like"] else "-u",
        OUTDIR_HAP="{sample_dir}/{sample_name}/{date}/references/haplotypes",
        SIM_BENCH=config.applications["simBench"],
    log:
        outfile="{sample_dir}/{sample_name}/{date}/references/haplotypes/simulate_haplotypes.out.log",
        errfile="{sample_dir}/{sample_name}/{date}/references/haplotypes/simulate_haplotypes.out.log",
    conda:
        config.simulate_haplotypes["conda"]
    shell:
        """
        if [[ -f "{params.HAPLOTYPE_SEQS}" ]]; then
            # use haplotypes from input sequence
            {params.SIM_BENCH} -f {params.HAPLOTYPE_SEQS} -n {params.NUM_HAPLOTYPES} -mr {params.MUT_RATE} -dr {params.DEL_RATE} -ir {params.INS_RATE} {params.NO_FR} {params.DEL_LEN} -s {params.SEED} -v -oh {params.OUTDIR_HAP} -o haplotypes > >(tee {log.outfile}) 2>&1
        else
            # use master sequence
            {params.SIM_BENCH} -f {input} {params.TREE_LIKE} -n {params.NUM_HAPLOTYPES} -mr {params.MUT_RATE} -dr {params.DEL_RATE} -ir {params.INS_RATE} {params.NO_FR} {params.DEL_LEN} -s {params.SEED} -v -oh {params.OUTDIR_HAP} -o haplotypes > >(tee {log.outfile}) 2>&1
        fi
        """


# 2. Simulate reads

if config.input["paired"]:

    rule simulate_reads:
        input:
            "{sample_dir}/{sample_name}/{date}/references/haplotypes/haplotypes.fasta",
        output:
            R1_raw="{sample_dir}/{sample_name}/{date}/raw_data/simreads_R1.fastq",
            R2_raw="{sample_dir}/{sample_name}/{date}/raw_data/simreads_R2.fastq",
        params:
            scratch="2000",
            mem=config.simulate_reads["mem"],
            time=config.simulate_reads["time"],
            NUM_HAPLOTYPES=lambda wildcards: sample_dict[
                sample_record(sample_name=wildcards.sample_name, date=wildcards.date)
            ]["num_haplotypes"],
            COVERAGE=lambda wildcards: sample_dict[
                sample_record(sample_name=wildcards.sample_name, date=wildcards.date)
            ]["coverage"],
            NUM_READS="-nr" if config.simulate_reads["num_reads"] else "",
            HIGH_QUAL="-q" if config.simulate_reads["high_quality"] else "",
            READ_LEN=lambda wildcards: sample_dict[
                sample_record(sample_name=wildcards.sample_name, date=wildcards.date)
            ]["read_len"],
            PAIRED="-p",
            FRAGMENT_SIZE=lambda wildcards: sample_dict[
                sample_record(sample_name=wildcards.sample_name, date=wildcards.date)
            ]["fragment_size"],
            FREQ_PARAMS=get_freq_params,
            FREQ_DSTR=lambda wildcards: sample_dict[
                sample_record(sample_name=wildcards.sample_name, date=wildcards.date)
            ]["freq_dstr"],
            SEED=lambda wildcards: sample_dict[
                sample_record(sample_name=wildcards.sample_name, date=wildcards.date)
            ]["seed"],
            OUTDIR_HAP="{sample_dir}/{sample_name}/{date}/references/haplotypes",
            OUTDIR_READS="{sample_dir}/{sample_name}/{date}/raw_data",
            ART=config.applications["art"],
            SIM_BENCH=config.applications["simBench"],
        log:
            outfile="{sample_dir}/{sample_name}/{date}/raw_data/simBench.out.log",
            errfile="{sample_dir}/{sample_name}/{date}/raw_data/simBench.out.log",
        conda:
            config.simulate_reads["conda"]
        shell:
            """
            {params.SIM_BENCH} -n {params.NUM_HAPLOTYPES} -c {params.COVERAGE} {params.NUM_READS} -l {params.READ_LEN} {params.PAIRED} -m {params.FRAGMENT_SIZE} -d {params.FREQ_DSTR} {params.FREQ_PARAMS} {params.HIGH_QUAL} -art {params.ART} -s {params.SEED} -v -oh {params.OUTDIR_HAP} -or {params.OUTDIR_READS} -o reads > >(tee {log.outfile}) 2>&1

            # Move intermediate results
            mkdir -p {params.OUTDIR_READS}/reads
            mv {params.OUTDIR_READS}/*.sam {params.OUTDIR_READS}/reads
            mv {params.OUTDIR_READS}/*.aln {params.OUTDIR_READS}/reads
            """


else:

    rule simulate_reads:
        input:
            "{sample_dir}/{sample_name}/{date}/references/haplotypes/haplotypes.fasta",
        output:
            R1_raw="{sample_dir}/{sample_name}/{date}/raw_data/simreads_R1.fastq",
        params:
            scratch="2000",
            mem=config.simulate_reads["mem"],
            time=config.simulate_reads["time"],
            NUM_HAPLOTYPES=lambda wildcards: sample_dict[
                sample_record(sample_name=wildcards.sample_name, date=wildcards.date)
            ]["num_haplotypes"],
            COVERAGE=lambda wildcards: sample_dict[
                sample_record(sample_name=wildcards.sample_name, date=wildcards.date)
            ]["coverage"],
            NUM_READS="-nr" if config.simulate_reads["num_reads"] else "",
            HIGH_QUAL="-q" if config.simulate_reads["high_quality"] else "",
            READ_LEN=lambda wildcards: sample_dict[
                sample_record(sample_name=wildcards.sample_name, date=wildcards.date)
            ]["read_len"],
            FRAGMENT_SIZE=lambda wildcards: sample_dict[
                sample_record(sample_name=wildcards.sample_name, date=wildcards.date)
            ]["fragment_size"],
            FREQ_PARAMS=get_freq_params,
            FREQ_DSTR=lambda wildcards: sample_dict[
                sample_record(sample_name=wildcards.sample_name, date=wildcards.date)
            ]["freq_dstr"],
            SEED=lambda wildcards: sample_dict[
                sample_record(sample_name=wildcards.sample_name, date=wildcards.date)
            ]["seed"],
            OUTDIR_HAP="{sample_dir}/{sample_name}/{date}/references/haplotypes",
            OUTDIR_READS="{sample_dir}/{sample_name}/{date}/raw_data",
            ART=config.applications["art"],
            SIM_BENCH=config.applications["simBench"],
        log:
            outfile="{sample_dir}/{sample_name}/{date}/raw_data/simBench.out.log",
            errfile="{sample_dir}/{sample_name}/{date}/raw_data/simBench.out.log",
        conda:
            config.simulate_reads["conda"]
        shell:
            """
            {params.SIM_BENCH} -n {params.NUM_HAPLOTYPES} -c {params.COVERAGE} {params.NUM_READS} -l {params.READ_LEN} -m {params.FRAGMENT_SIZE} -d {params.FREQ_DSTR} {params.FREQ_PARAMS} {params.HIGH_QUAL} -art {params.ART} -s {params.SEED} -v -oh {params.OUTDIR_HAP} -or {params.OUTDIR_READS} -o reads > >(tee {log.outfile}) 2>&1

            # Move intermediate results
            mkdir -p {params.OUTDIR_READS}/reads
            mv {params.OUTDIR_READS}/*.sam {params.OUTDIR_READS}/reads
            mv {params.OUTDIR_READS}/*.aln {params.OUTDIR_READS}/reads
            """
