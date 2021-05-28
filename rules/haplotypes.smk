__author__ = "Susana Posada-Cespedes"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"


rule haploclique:
    input:
        "{dataset}/alignments/REF_aln.bam",
    output:
        FASTA="{dataset}/variants/global/quasispecies.fasta",
        BAM="{dataset}/variants/global/quasispecies.bam",
    params:
        scratch="1250",
        mem=config.haploclique["mem"],
        time=config.haploclique["time"],
        RELAX=(
            "--edge_quasi_cutoff_cliques=0.85 --edge_quasi_cutoff_mixed=0.85 --edge_quasi_cutoff_single=0.8 --min_overlap_cliques=0.6 --min_overlap_single=0.5"
            if config.haploclique["relax"]
            else ""
        ),
        NO_SINGLETONS="--no_singletons" if config.haploclique["no_singletons"] else "",
        NO_PROB0="--no_prob0" if config.haploclique["no_prob0"] else "",
        CLIQUE_SIZE_LIMIT=config.haploclique["clique_size_limit"],
        MAX_NUM_CLIQUES=config.haploclique["max_num_cliques"],
        OUTPREFIX="{dataset}/variants/global/quasispecies",
        HAPLOCLIQUE=config.applications["haploclique"],
    log:
        outfile="{dataset}/variants/global/haploclique.out.log",
        errfile="{dataset}/variants/global/haploclique.err.log",
    conda:
        config.haploclique["conda"]
    benchmark:
        "{dataset}/variants/global/haploclique.benchmark"
    threads: 1
    shell:
        """
        {params.HAPLOCLIQUE} {params.RELAX} {params.NO_SINGLETONS} {params.NO_PROB0} --limit_clique_size={params.CLIQUE_SIZE_LIMIT} --max_cliques={params.MAX_NUM_CLIQUES} --log={log.outfile} --bam {input} {params.OUTPREFIX} 2> >(tee {log.errfile} >&2)
        """


rule haploclique_visualization:
    input:
        BAM="{dataset}/variants/global/quasispecies.bam",
        FASTA="{dataset}/variants/global/quasispecies.fasta",
    output:
        PDF="{dataset}/variants/global/quasispecies_plot.pdf",
    params:
        scratch="1250",
        mem=config.haploclique_visualization["mem"],
        time=config.haploclique_visualization["time"],
        REGION_START=config.haploclique_visualization["region_start"],
        REGION_END=config.haploclique_visualization["region_end"],
        USE_MSA="-r" if len(config.haploclique_visualization["msa"]) > 0 else "",
        MSA=config.haploclique_visualization["msa"],
        TSV="{dataset}/variants/global/quasispecies_mapping.tsv",
        INPREFIX="{dataset}/variants/global/quasispecies",
        COMPUTE_MDS=config.applications["compute_mds"],
    log:
        outfile="{dataset}/variants/global/haploclique_visualization.out.log",
        errfile="{dataset}/variants/global/haploclique_visualization.err.log",
    conda:
        config.haploclique_visualization["conda"]
    benchmark:
        "{dataset}/variants/global/haploclique_visualization.benchmark"
    threads: 1
    shell:
        """
        {params.COMPUTE_MDS} -q {params.INPREFIX} -s {params.REGION_START} -e {params.REGION_END} {params.USE_MSA} {params.MSA} -p {output.PDF} -o {params.TSV} > {log.output} 2> >(tee {log.errfile} >&2)
        """


if config.input["paired"]:

    rule savage:
        input:
            "{dataset}/alignments/REF_aln.bam",
        output:
            R1=temp("{dataset}/variants/global/R1.fastq"),
            R2=temp("{dataset}/variants/global/R2.fastq"),
            FASTA="{dataset}/variants/global/contigs_stage_c.fasta",
        params:
            scratch="1250",
            mem=config.savage["mem"],
            time=config.savage["time"],
            SPLIT=config.savage["split"],
            PICARD=config.applications["picard"],
            SAVAGE=config.applications["savage"],
            OUTDIR="{dataset}/variants/global/",
            FUNCTIONS=functions,
        log:
            outfile="{dataset}/variants/global/savage.out.log",
            errfile="{dataset}/variants/global/savage.err.log",
        conda:
            config.savage["conda"]
        benchmark:
            "{dataset}/variants/global/savage.benchmark"
        threads: config.savage["threads"]
        shell:
            """
            # Convert BAM to FASTQ without re-reversing reads - SAVAGE expect all reads in the same direction
            source {params.FUNCTIONS}
            SamToFastq {params.PICARD} I={input} FASTQ={output.R1} SECOND_END_FASTQ={output.R2} RC=false 2> >(tee -a {log.errfile} >&2)
            # Remove /1 and /2 from the read names
            sed -i -e "s:/1$::" {output.R1}
            sed -i -e "s:/2$::" {output.R2}

            R1=${{PWD}}/{output.R1}
            R2=${{PWD}}/{output.R2}
            {params.SAVAGE} -t {threads} --split {params.SPLIT} -p1 ${{R1}} -p2 ${{R2}} -o {params.OUTDIR} 2> >(tee -a {log.errfile} >&2)
            """


else:

    rule savage_se:
        input:
            "{dataset}/alignments/REF_aln.bam",
        output:
            R1=temp("{dataset}/variants/global/R1.fastq"),
            FASTA="{dataset}/variants/global/contigs_stage_c.fasta",
        params:
            scratch="1250",
            mem=config.savage["mem"],
            time=config.savage["time"],
            SPLIT=config.savage["split"],
            PICARD=config.applications["picard"],
            SAVAGE=config.applications["savage"],
            OUTDIR="{dataset}/variants/global/",
            FUNCTIONS=functions,
        log:
            outfile="{dataset}/variants/global/savage.out.log",
            errfile="{dataset}/variants/global/savage.err.log",
        conda:
            config.savage["conda"]
        benchmark:
            "{dataset}/variants/global/savage.benchmark"
        threads: config.savage["threads"]
        shell:
            """
            # Convert BAM to FASTQ without re-reversing reads - SAVAGE expect all reads in the same direction
            source {params.FUNCTIONS}
            SamToFastq {params.PICARD} I={input} FASTQ={output.R1} 2> >(tee -a {log.errfile} >&2)

            R1=${{PWD}}/{output.R1}
            {params.SAVAGE} -t {threads} --split {params.SPLIT} -s ${{R1}} -o {params.OUTDIR} 2> >(tee -a {log.errfile} >&2)
            """
