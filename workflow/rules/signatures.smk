import os

__author__ = "Ivan Topolsky"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"


def all_vocs(dir):
    (VOCS,) = glob_wildcards(os.path.join(dir, "{voc}.yaml"))
    return expand(os.path.join(dir, "{voc}.yaml"), voc=VOCS)


def proto_inserts(wildcards):
    return protocol_proto_option(proto=wildcards.proto, option="inserts_bedfile")


rule amplicons:
    input:
        vocs=all_vocs("references/voc/"),
        inserts=proto_inserts,
    output:
        amplicons=cohortdir("amplicons.{proto}.yaml"),
    params:
        COJAC="../cojac/cojac-wrapper",
        vocdir="references/voc/",
    log:
        outfile=cohortdir("amplicons.{proto}.out.log"),
        errfile=cohortdir("amplicons.{proto}.err.log"),
    conda:
        "../envs/cojac.yaml"
    benchmark:
        cohortdir("amplicons.{proto}.benchmark")
    resources:
        disk_mb=1024,
        mem_mb=256,
        runtime=10,
    threads: 1
    shell:
        """
        vocs=( {input.vocs} )
        {params.COJAC} cooc-mutbamscan "${{vocs[@]/#/--voc=}}" --bedfile="{input.inserts}" --out-amplicons="{output.amplicons}"  2> >(tee -a {log.errfile} >&2)  > >(tee -a {log.outfile})
        """


def get_sample_amplicons(wildcards):
    return cohortdir("amplicons.%s.yaml" % get_sample_protocol(wildcards))


rule cooc:
    input:
        BAM=alignment_wildcard,
        amplicons=get_sample_amplicons,
    output:
        cooc_yaml="{dataset}/signatures/cooc.yaml",
        cooc_csv="{dataset}/signatures/cooc.csv",
    params:
        COJAC="../cojac/cojac-wrapper",
        name=ID,
        sep=config.general["id_separator"],
    log:
        outfile="{dataset}/signatures/cooc.out.log",
        errfile="{dataset}/signatures/cooc.err.log",
    conda:
        "../envs/cojac.yaml"
    benchmark:
        "{dataset}/signatures/cooc.benchmark"
    resources:
        disk_mb=1024,
        mem_mb=8192,
        runtime=45,
    threads: 1
    shell:
        """
        {params.COJAC} cooc-mutbamscan --alignments="{input.BAM}" --name="{params.name}" --in-amp="{input.amplicons}" --yaml="{output.cooc_yaml}"   2> >(tee -a {log.errfile} >&2)  > >(tee -a {log.outfile})
        {params.COJAC} cooc-tabmut --yaml="{output.cooc_yaml}" --output="{output.cooc_csv}" --multiindex --lines --batchname="{params.sep}" 2> >(tee -a {log.errfile} >&2)  > >(tee -a {log.outfile})
        """


def get_s_rec(wildcards):
    return guess_sample(wildcards.dataset)


rule sigmut:
    input:
        basecnt="{dataset}/alignments/basecnt.tsv.gz",
        mutlist="mutlist.txt",
    output:
        mut="{dataset}/signatures/mut.tsv.gz",
    params:
        getmutation=cachepath(
            "../scripts/getmutation_from_sample.py",
            executable=True,
            localsource=True,
        ),
        s_rec=get_s_rec,
    log:
        outfile="{dataset}/signatures/mut.out.log",
        errfile="{dataset}/signatures/mut.err.log",
    conda:
        "../envs/cojac.yaml"
    benchmark:
        "{dataset}/signatures/mut.benchmark"
    resources:
        disk_mb=1024,
        mem_mb=1024,
        runtime=15,
    threads: 1
    shell:
        """
        {params.getmutation} --outname {output.mut} --sample {params.s_rec.sample_id} --batch {params.s_rec.date} -m {input.mutlist} -- {input.basecnt} 2> >(tee -a {log.errfile} >&2)  > >(tee -a {log.outfile})
        """


rule timeline:
    input:
        "samples.wastewateronly.tsv",
    output:
        "timeline.tsv",
    params:
        locations="ww_locations.tsv",
        regex="regex.yaml",
        maketimeline=cachepath(
            "../scripts/file_parser.py",
            executable=True,
            localsource=True,
        ),
    log:
        outfile=cohortdir("timeline.out.log"),
        errfile=cohortdir("timeline.err.log"),
    conda:
        "../envs/lollipop.yaml"
    benchmark:
        cohortdir("timeline.err.benchmark")
    resources:
        disk_mb=1024,
        mem_mb=1024,
        runtime=15,
    threads: 1
    shell:
        """
        {params.maketimeline} --regex-config {params.regex} --no-fallback --locations {params.locations} --output {output} -- {input} 2> >(tee -a {log.errfile} >&2)  > >(tee -a {log.outfile})
        """
