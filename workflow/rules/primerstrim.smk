from collections import UserDict
import re

__author__ = "Michal Okoniewski"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"


# Assuming that getting the indexed BAMs from sam2bam in align.smk (ca l 427)


# fetch primers_file either from config or sample-specific protocol
re_trim_dataset = re.compile("(?P<dataset>.+)/alignments/REF_aln")


def primers_file(wildcards):
    # skip if no sample ever has 4th column
    if 0 == sample_proto_count:
        return config.input["primers_bedfile"]

    # file -> dataset
    ds = re_trim_dataset.match(wildcards.file)
    if not ds:
        raise ValueError(f"cannot guess the sample dataset for file {wildcards.file}")

    wc = UserDict(ds.groupdict())
    wc.__dict__.update(wc)
    return protocol_option(wc, option="primers_bedfile")


rule primerstrim:
    input:
        BAM="{file}.bam",
        BAI="{file}.bam.bai",
    output:
        BAM="{file}_trim.bam",
        BAI="{file}_trim.bam.bai",
    params:
        SAMTOOLS=config.applications["samtools"],
        IVAR=config.applications["ivar"],
        BED_PRIMERS=primers_file,
        # prefixes to use for both softwares
        ivar_tmp=temp_prefix("{file}_trim_unsorted"),
        sort_tmp=temp_prefix("{file}_trim_tmp"),
    log:
        outfile="{file}_trim.out.log",
        errfile="{file}_trim.err.log",
    conda:
        config.primerstrim["conda"]
    benchmark:
        "{file}_trim.benchmark"
    resources:
        disk_mb=1250,
        mem_mb=config.primerstrim["mem"],
        runtime=config.primerstrim["time"],
    threads: 1
    shell:
        """
        echo "Trimming BAM with ivar"

        # iVar will Segfault without this:
        mkdir -p "$(dirname {params.ivar_tmp}"")"
        {params.IVAR} trim -e -i {input.BAM} -b {params.BED_PRIMERS} -p {params.ivar_tmp} > {log.outfile} 2> {log.errfile}

        # samtools complains without that:
        rm -f '{params.sort_tmp}'.[0-9]*.bam
        {params.SAMTOOLS}  sort -o {output.BAM} -T {params.sort_tmp} {params.ivar_tmp}.bam 2> >(tee -a {log.errfile} >&2)
        {params.SAMTOOLS}  index {output.BAM} 2> >(tee -a {log.errfile} >&2)
        rm -f {params.ivar_tmp}.bam

        """


rule ampliconclip:
    input:
        BAM="{file}.bam",
        BAI="{file}.bam.bai",
    output:
        BAM="{file}_trim.bam",
        BAI="{file}_trim.bam.bai",
        stats="{file}_trim.stats",
    params:
        SAMTOOLS=config.applications["samtools"],
        BED_PRIMERS=primers_file,
        # prefixes to use for both softwares
        sort_tmp=temp_prefix("{file}_trim_tmp"),
    log:
        outfile="{file}_trim.out.log",
        errfile="{file}_trim.err.log",
    conda:
        config.sam2bam["conda"]
    benchmark:
        "{file}_trim.benchmark"
    resources:
        disk_mb=1250,
        mem_mb=config.primerstrim["mem"],
        runtime=config.primerstrim["time"],
    threads: 1
    shell:
        """
        echo "Trimming BAM with samtools"

        # samtools complains without that:
        rm -f '{params.sort_tmp}'.[0-9]*.bam

        {params.SAMTOOLS} ampliconclip -b {params.BED_PRIMERS} -f {output.stats} {input.BAM} |
            {params.SAMTOOLS}  sort -o {output.BAM} -T {params.sort_tmp} 2> >(tee {log.errfile} >&2)
        {params.SAMTOOLS}  index {output.BAM} 2> >(tee -a {log.errfile} >&2)
        """


if config.general["primers_trimmer"] == "ivar":

    ruleorder: primerstrim > ampliconclip


elif config.general["primers_trimmer"] == "samtools":

    ruleorder: ampliconclip > primerstrim
