![Logo](https://cbg-ethz.github.io/V-pipe/img/logo.svg)

V-pipe is a workflow designed for analysis of next generation sequencing (NGS) data from viral pathogens. It produces a number of results in a curated format.

[![Snakemake](https://img.shields.io/badge/snakemake-≥4.8.0-brightgreen.svg?style=plastic)](https://snakemake.bitbucket.io)
[![License: Apache-2.0](https://img.shields.io/badge/License-Apache_2.0-yellow.svg?style=plastic)](https://opensource.org/licenses/Apache-2.0)
[![Analytics](https://cgb-ethz-ga-beacon.appspot.com/UA-126950686-1/github/README)](https://github.com/igrigorik/ga-beacon)

## Dependencies

- **Snakemake**

  Snakemake is the central workflow and dependency manager of V-pipe. It determines the order in which individual tools are invoked and checks that programs do not exit unexpectedly.

- **FastQC**

  FastQC gives an overview of the raw sequencing data. Flowcells that have been overloaded or otherwise fail during sequencing can easily be determined with FastQC.

- **PRINSEQ**

  Trimming and clipping of reads is performed by PRINSEQ. It is currently the most versatile raw read processor with many customization options.

- **Vicuna**

  Vicuna is a de novo assembler designed for generating rough reference contigs of viral NGS data. It can deal with the inherent heterogeneity such as high single-base heterogeneity and structural variants.

- **ngshmmalign**

  We perform the alignment of the curated NGS data using our custom ngshmmalign that takes structural variants into account. It produces multiple consensus sequences that include either majority bases or ambiguous bases.

- **bwa**

  In order to detect specific cross-contaminations with other probes, the Burrows-Wheeler aligner is used. It quickly yields estimates for foreign genomic material in an experiment.

- **MAFFT**

  To standardise multiple samples to the same reference genome (say HXB2 for HIV-1), the multiple sequence aligner MAFFT is employed. The multiple sequence alignment helps in determining regions of low conservation and thus makes standardisation of alignments mroe robust.

- **QuasiRecomb**

  QuasiRecomb performs local and global haplotype reconstruction for heterogeneous viral populations by using a hidden Markov model.

- **Samtools**

  The Swiss Army knife of alignment postprocessing and diagnostics.

- **SmallGenomeUtilities**

  We perform genomic liftovers to standardised reference genomes using our in-house developed python library of utilities for rewriting alignments.

## Contributions

- [David Seifert](https://orcid.org/0000-0003-4739-5110)
- [Susana Posada Céspedes](https://orcid.org/0000-0002-7459-8186)
- [Niko Beerenwinkel](https://orcid.org/0000-0002-0573-6119)

## Contacts

For any request regarding this pipeline, contact the V-pipe Dev Team <v-pipe@bsse.ethz.ch>
