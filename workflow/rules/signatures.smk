import os

__author__ = "Ivan Topolsky"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"


def all_vocs(dir):
    if not dir:
        return []
    (VOCS,) = glob_wildcards(os.path.join(dir, "{voc}.yaml"))
    return expand(os.path.join(dir, "{voc}.yaml"), voc=VOCS)


def proto_inserts(wildcards):
    return protocol_proto_option(proto=wildcards.proto, option="inserts_bedfile")


rule amplicons:
    input:
        vocs=all_vocs(config.input["variants_def_directory"]),
        inserts=proto_inserts,
    output:
        amplicons=cohortdir("amplicons.{proto}.yaml"),
    params:
        COJAC=config.applications["cojac"],
        vocdir=config.input["variants_def_directory"],
        mincooc=config.amplicons["mincooc"],
    log:
        outfile=cohortdir("amplicons.{proto}.out.log"),
        errfile=cohortdir("amplicons.{proto}.err.log"),
    conda:
        config.cooc["conda"]
    benchmark:
        cohortdir("amplicons.{proto}.benchmark")
    resources:
        disk_mb=1024,
        mem_mb=config.amplicons["mem"],
        runtime=config.amplicons["time"],
    threads: 1
    shell:
        """
        vocs=( {input.vocs} )
        {params.COJAC} cooc-mutbamscan "${{vocs[@]/#/--voc=}}" --bedfile="{input.inserts}" --cooc="{params.mincooc}" --out-amplicons="{output.amplicons}"  2> >(tee -a {log.errfile} >&2)  > >(tee -a {log.outfile})
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
        COJAC=config.applications["cojac"],
        name=ID,
        sep=config.general["id_separator"],
        out_format=(
            "--multiindex"
            if config.cooc["out_format"] == "columns"
            else "--multiindex --lines"
        ),
    log:
        outfile="{dataset}/signatures/cooc.out.log",
        errfile="{dataset}/signatures/cooc.err.log",
    conda:
        config.cooc["conda"]
    benchmark:
        "{dataset}/signatures/cooc.benchmark"
    resources:
        disk_mb=1024,
        mem_mb=config.cooc["mem"],
        runtime=config.cooc["time"],
    threads: config.cooc["threads"]
    shell:
        """
        {params.COJAC} cooc-mutbamscan --alignments="{input.BAM}" --name="{params.name}" --in-amp="{input.amplicons}" --yaml="{output.cooc_yaml}"   2> >(tee -a {log.errfile} >&2)  > >(tee -a {log.outfile})
        {params.COJAC} cooc-tabmut --yaml="{output.cooc_yaml}" --output="{output.cooc_csv}" {params.out_format} --batchname="{params.sep}" 2> >(tee -a {log.errfile} >&2)  > >(tee -a {log.outfile})
        """


def proto_datasets(wildcards):
    return [
        sample_paths[s_rec]
        for s_rec, s_row in sample_table.items()
        if s_row.protocol == wildcards.proto
    ]


def cohort_cooc_proto(wildcards):
    return expand("{dataset}/signatures/cooc.yaml", dataset=proto_datasets(wildcards))


rule cohort_cooc:
    input:
        cohort_cooc_proto,
    output:
        cooc_yaml=cohortdir("cohort_cooc.{proto}.yaml"),
        cooc_csv=cohortdir("cohort_cooc.{proto}.csv"),
    params:
        COJAC=config.applications["cojac"],
        sep=config.general["id_separator"],
        out_format=(
            "--multiindex"
            if config.cooc["out_format"] == "columns"
            else "--multiindex --lines"
        ),
    log:
        outfile=cohortdir("cohort_cooc.{proto}.out.log"),
        errfile=cohortdir("cohort_cooc.{proto}.err.log"),
    conda:
        config.cooc["conda"]
    benchmark:
        cohortdir("cohort_cooc.{proto}.benchmark")
    resources:
        disk_mb=1024,
        mem_mb=config.cooc["mem"],
        runtime=config.cooc["time"],
    threads: config.cooc["threads"]
    shell:
        """
        cat {input} > {output.cooc_yaml} 2> >(tee {log.errfile} >&2)
        {params.COJAC} cooc-tabmut --yaml="{output.cooc_yaml}" --output="{output.cooc_csv}" {params.out_format} --batchname="{params.sep}" 2> >(tee -a {log.errfile} >&2)  > >(tee {log.outfile})
        """


rule cohort_cooc_report:
    input:
        cooc_yaml=cohortdir("cohort_cooc.{proto}.yaml"),
        amplicons=cohortdir("amplicons.{proto}.yaml"),
    output:
        cooc_report_csv=cohortdir("cohort_cooc_report.{proto}.csv"),
    params:
        COJAC=config.applications["cojac"],
        vocdir=config.input["variants_def_directory"],
        sep=config.general["id_separator"],
    log:
        outfile=cohortdir("cohort_cooc_report.{proto}.out.log"),
        errfile=cohortdir("cohort_cooc_report.{proto}.err.log"),
    conda:
        config.cooc["conda"]
    benchmark:
        cohortdir("cohort_cooc_report.{proto}.benchmark")
    resources:
        disk_mb=1024,
        mem_mb=config.cooc["mem"],
        runtime=config.cooc["time"],
    threads: config.cooc["threads"]
    shell:
        """
        {params.COJAC} cooc-pubmut --yaml="{input.cooc_yaml}" --amplicons="{input.amplicons}" --vocdir="{params.vocdir}" --output="{output.cooc_report_csv}" --batchname="{params.sep}" 2> >(tee {log.errfile} >&2)  > >(tee {log.outfile})
        """


rule mutlist:
    input:
        vocs=all_vocs(config.input["variants_def_directory"]),
        gff=config.input["genes_gff"],
    output:
        mutlist=cohortdir("mutlist.tsv"),
        pagovars=cohortdir("variants_pangolin.yaml"),
    params:
        LOLLIPOP=config.applications["lollipop"],
    log:
        outfile=cohortdir("mutlist.out.log"),
        errfile=cohortdir("mutlist.err.log"),
    conda:
        config.deconvolution["conda"]
    benchmark:
        cohortdir("mutlist.benchmark")
    resources:
        disk_mb=1024,
        mem_mb=config.mutlist["mem"],
        runtime=config.mutlist["time"],
    threads: 1
    shell:
        """
        {params.LOLLIPOP} generate-mutlist --output {output.mutlist} --out-pangovars {output.pagovars} --genes {input.gff} -- {input.vocs}
        """


def get_s_rec(wildcards):
    return guess_sample(wildcards.dataset)


rule sigmut:
    input:
        basecnt="{dataset}/alignments/basecnt.tsv.gz",
        mutlist=cohortdir("mutlist.tsv"),
    output:
        mut="{dataset}/signatures/mut.tsv",
    params:
        LOLLIPOP=config.applications["lollipop"],
        ARRAYBASED=config.general["tsvbased"],
        s_rec=get_s_rec,
    log:
        outfile="{dataset}/signatures/mut.out.log",
        errfile="{dataset}/signatures/mut.err.log",
    conda:
        config.deconvolution["conda"]
    benchmark:
        "{dataset}/signatures/mut.benchmark"
    resources:
        disk_mb=1024,
        mem_mb=config.sigmut["mem"],
        runtime=config.sigmut["time"],
    threads: 1
    shell:
        """
        {params.LOLLIPOP} getmutations from-basecount --outname "{output.mut}" --samplename "{params.s_rec.sample_id}" --batch "{params.s_rec.date}" -m "{input.mutlist}" --based "{params.ARRAYBASED}" -- "{input.basecnt}" 2> >(tee -a {log.errfile} >&2)  > >(tee -a {log.outfile})
        """


if config.timeline["local"]:

    localrules:
        timeline,


rule timeline:
    input:
        samples_tsv=config.input["samples_file"],
        locations=(
            config.timeline["locations_table"]
            if config.timeline["locations_table"]
            else []
        ),
        regex=config.timeline["regex_yaml"] if config.timeline["regex_yaml"] else [],
    output:
        timeline=cohortdir("timeline.tsv"),
        locations_list=(
            cohortdir("locations_list.yaml")
            if config.timeline["locations_table"]
            else []
        ),
    params:
        maketimeline=cachepath(config.timeline["script"], executable=True),
        locations=(
            f"--locations {config.timeline['locations_table']}"
            if config.timeline["locations_table"]
            else ""
        ),
        regex=(
            f"--regex-config {config.timeline['regex_yaml']}"
            if config.timeline["regex_yaml"]
            else ""
        ),
        out_locations=(
            f"--out-locations {cohortdir('locations_list.yaml')}"
            if config.timeline["locations_table"] or config.timeline["regex_yaml"]
            else ""
        ),
        options=config.timeline["options"],
    log:
        outfile=cohortdir("timeline.out.log"),
        errfile=cohortdir("timeline.err.log"),
    conda:
        config.timeline["conda"]
    benchmark:
        cohortdir("timeline.benchmark")
    resources:
        disk_mb=1024,
        mem_mb=config.timeline["mem"],
        runtime=config.timeline["time"],
    threads: config.timeline["threads"]
    shell:
        """
        {params.maketimeline} {params.regex} {params.locations} {params.out_locations} --output "{output.timeline}" {params.options} -- "{input.samples_tsv}" 2> >(tee -a {log.errfile} >&2)  > >(tee -a {log.outfile})
        """


rule tallymut:
    input:
        muts=expand("{dataset}/signatures/mut.tsv", dataset=datasets),
        times=(config.tallymut.get("timeline_file", None) or cohortdir("timeline.tsv")),
    output:
        tallymut=cohortdir("tallymut.tsv.zst"),
    params:
        XSV=config.applications["xsv"],
        ZSTD=config.applications["zstd"],
        selector="sample,batch" if sample_2level_count else "sample",
    log:
        outfile=cohortdir("tallymut.out.log"),
        errfile=cohortdir("tallymut.err.log"),
    conda:
        config.tallymut["conda"]
    benchmark:
        cohortdir("tallymut.benchmark")
    resources:
        disk_mb=1024,
        mem_mb=config.tallymut["mem"],
        runtime=config.tallymut["time"],
    threads: 1
    shell:
        """
        {params.XSV} join --right {params.selector} {input.times} {params.selector} \
         <({params.XSV} cat rows --delimiter '\\t' {input.muts} ) \
         | {params.XSV} fmt --out-delimiter '\\t' \
         | {params.ZSTD} -o {output.tallymut} 2> >(tee -a {log.errfile} >&2) > >(tee -a {log.outfile})
        """


rule deconvolution:
    input:
        tallymut=cohortdir("tallymut.tsv.zst"),
        deconv_conf=config.deconvolution["deconvolution_config"],
        var_conf=(
            config.deconvolution["variants_config"]
            if config.deconvolution["variants_config"]
            else cohortdir("variants_pangolin.yaml")
        ),
        var_dates=(
            config.deconvolution["variants_dates"]
            if config.deconvolution["variants_dates"]
            else []
        ),
        filters=(
            config.deconvolution["filters"] if config.deconvolution["filters"] else []
        ),
    output:
        deconvoluted=cohortdir("deconvoluted.tsv.zst"),
        deconv_json=cohortdir("deconvoluted_upload.json"),
    params:
        LOLLIPOP=config.applications["lollipop"],
        out_format=(
            "--fmt-columns" if config.deconvolution["out_format"] == "columns" else ""
        ),
        seed="--seed=42",
    log:
        outfile=cohortdir("deconvoluted.out.log"),
        errfile=cohortdir("deconvoluted.err.log"),
    conda:
        config.deconvolution["conda"]
    benchmark:
        cohortdir("deconvoluted.benchmark")
    resources:
        disk_mb=1024,
        mem_mb=config.deconvolution["mem"],
        runtime=config.deconvolution["time"],
    threads: config.deconvolution["threads"]
    shell:
        """
        {params.LOLLIPOP} deconvolute "--output={output.deconvoluted}" "--out-json={output.deconv_json}" "--var={input.var_conf}" "--vd={input.var_dates}" "--dec={input.deconv_conf}" "--filters={input.filters}" {params.out_format} {params.seed} "{input.tallymut}" 2> >(tee -a {log.errfile} >&2) > >(tee -a {log.outfile})
        """


rule allCooc:
    input:
        expand("{dataset}/signatures/cooc.yaml", dataset=datasets),
