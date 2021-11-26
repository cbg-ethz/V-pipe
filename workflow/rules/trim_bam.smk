import os

__author__ = "Michal Okoniewski"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"


os.system(shellcommand)
shellcommand="mkdir vtrim"

# Assuming that getting the indexed BAMs from sam2bam in align.smk (ca l 427)


rule ivar_trim:
    input:
        BAM="{file}.bam",
        BAI="{file}.bam.bai",
    output:
        BAM="vtrim/{file}.bam",
        BAI="vtrim/{file}.bam.bai"
    params:
        SAMTOOLS=config.applications["samtools"],
        FUNCTIONS=functions,
        BED_PRIMERS="resources/sars-cov-2/primers/nCoV-2019.primer.bed",
    conda:
       config.ivar["conda"],
    threads: 1
    shell:
        """
        echo "Trimming BAM with ivar"

        ivar trim -e -i {input.BAM} -b {params.BED_PRIMERS} -p vtrim/unsorted_{wildcards.file}
        {params.SAMTOOLS}  sort -o {output.BAM} -T out/{wildcards.file} {input.BAM}
        {params.SAMTOOLS}  index {output.BAM}
        rm vtrim/unsorted_{wildcards.file}.bam

        """


