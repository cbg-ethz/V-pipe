# GROUP: local
# CONDA: lofreq = 2.1.5
# CONDA: pysam = 0.19.0


import subprocess
from pathlib import Path
import pysam

def make_quality_scores_uninformative(fname_bam_in, fname_bam_out):

    bamfile = pysam.AlignmentFile(fname_bam_in, "r")
    new_head = bamfile.header.to_dict()

    with pysam.AlignmentFile(fname_bam_out, "w", header=new_head) as outf:
        for read in bamfile.fetch():
            a = pysam.AlignedSegment(outf.header)
            a.query_name = read.query_name
            a.query_sequence = read.query_sequence
            a.reference_name = read.reference_name
            a.flag = read.flag
            a.reference_start = read.reference_start
            a.mapping_quality = read.mapping_quality
            a.cigar = read.cigar
            a.next_reference_id = read.next_reference_id
            a.next_reference_start = read.next_reference_start
            a.template_length = read.template_length
            quality_not_stored = ''.join(len(read.query_qualities)*['*'])
            a.query_qualities = pysam.qualitystring_to_array(quality_not_stored)
            a.tags = read.tags
            outf.write(a)



def main(fname_bam, fname_reference, fname_result, fname_result_haplo, dname_work):
    dname_work.mkdir(parents=True, exist_ok=True)


    fname_bam_in = fname_bam
    fname_bam_out = fname_bam.resolve().split('.bam')[0]+'NQ.bam'
    make_quality_scores_uninformative(fname_bam_in, fname_bam_out)

    subprocess.run(
        [
            "lofreq",
            "call",
            "-f",
            fname_reference.resolve(),
            "-o",
            fname_result.resolve(),
            fname_bam.resolve(),
        ],
        cwd=dname_work,
        check=True,
    )

    open(fname_result_haplo.resolve(), 'a').close()


if __name__ == "__main__":
    main(
        Path(snakemake.input.fname_bam),
        Path(snakemake.input.fname_reference),
        Path(snakemake.output.fname_result_snv),
        Path(snakemake.output.fname_result_haplo),
        Path(snakemake.output.dname_work),
    )
