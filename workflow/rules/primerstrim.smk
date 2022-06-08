__author__ = "Michal Okoniewski"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"


# Assuming that getting the indexed BAMs from sam2bam in align.smk (ca l 427)


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
        # TODO make it configurable
        BED_PRIMERS=config.input["primers_file"],
        # prefixes to use for both softwares
        ivar_tmp=temp_prefix("{file}_trim_unsorted"),
        sort_tmp=temp_prefix("{file}_trim_tmp"),
    log:
        outfile="{file}_trim.out.log",
        errfile="{file}_trim.err.log",
    conda:
        config.primerstrim["conda"]
    resources:
        disk_mb=1250,
        mem_mb=config.primerstrim["mem"],
        time_min=config.primerstrim["time"],
    threads: 1
    shell:
        """
        echo "Trimming BAM with ivar"

        {params.IVAR} trim -e -i {input.BAM} -b {params.BED_PRIMERS} -p {params.ivar_tmp}
        rm -f '{params.sort_tmp}'.[0-9]*.bam
        {params.SAMTOOLS}  sort -o {output.BAM} -T {params.sort_tmp} {params.ivar_tmp}.bam
        {params.SAMTOOLS}  index {output.BAM}
        rm -f {params.ivar_tmp}.bam

        """
