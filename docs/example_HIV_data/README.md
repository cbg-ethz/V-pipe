Sample HIV data for tutorials
=============================

The samples were taken from the publication Abrahams et al. (2019), Science translational medicine 11.513 (DOI: 10.1126/scitranslmed.aaw5589).

We download the following HIV Multiplexed Illumina MiSeq data from the short read archive (SRA): SRR9588830, SRR9588828, SRR9588844 and SRR9588785. They where taken from a HIV-1 positive patients at different time points post-infection.

To download the data, we used ``` sra-tools ```.

```bash
mkdir -p samples/CAP188/4/
cd samples/CAP188/4/
fastq-dump -O raw_data --split-e  SRR9588828
```

Using the `--split-e` option, we download the reads seperatled into forward and reverse reads. Here you can find some more information on `fastq-dump`: <https://edwards.flinders.edu.au/fastq-dump/>

We aligned the reads to the HIV strain HXB2, retrieved reads covering the region HXB2:2453-3356, and further subsampled to have a feasible sized sample to run on a laptop.

```bash
samtools view REF_aln.bam -h "HXB2:2453-3356" > output_region.bam
samtools view -s 0.10 -b output_region.bam > output_region_subsample.bam
```
