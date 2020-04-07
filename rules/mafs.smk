
# 1. Ouptut minor allele frequencies
rule minor_variants:
    input:
        REF = reference_file,
        BAM = expand("{dataset}/alignments/REF_aln.bam", dataset=datasets),
    output:
        VARIANTS = "variants/minority_variants.tsv",
        CONSENSUS = "variants/cohort_consensus.fasta",
        COVERAGE = "variants/coverage.tsv"
    params:
        scratch = '1250',
        mem = config.minor_variants['mem'],
        time = config.minor_variants['time'],
        OUTDIR = "variants",
        NAMES = IDs,
        MINORITY_CALLER = config.applications['minority_freq'],
    log:
        outfile = "variants/minority_variants.out.log",
        errfile = "variants/minority_variants.out.log",
    conda:
        config.minor_variants['conda']
    benchmark:
        "variants/minority_variants.benchmark"
    threads:
        config.minor_variants['threads']
    shell:
        """
        {params.MINORITY_CALLER} -r {input.REF} -N {params.NAMES} -t {threads} -o {params.OUTDIR} -d {input.BAM} > >(tee {log.outfile}) 2>&1
        """

