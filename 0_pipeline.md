---
layout: page
title: Pipeline overview
permalink: /pipeline/
---


![Pipeline](../img/pipeline.svg)


## Dependencies

- **[Conda](https://conda.io/docs/index.html)**
  [![Shell installer](https://img.shields.io/badge/bash-installer-lightgrey.svg?style=flat)](https://conda.io/miniconda.html)

  Conda is an open source package management system and environment management system. V-pipe uses it to automatically obtain reproducible environments and simplify installation of the individual components of the pipeline, thanks to the [Bioconda channel](https://bioconda.github.io) - a distribution of bioinformatics software.

  See the [documentation](http://conda.io/docs/install/quick.html) of conda to install it.

- **[Snakemake](https://snakemake.github.io/)**
  [![Bioconda package](https://img.shields.io/conda/dn/bioconda/snakemake.svg?label=Bioconda)](https://bioconda.github.io/recipes/snakemake/README.html)
  [![Snakemake](https://img.shields.io/badge/snakemake-≥4.8.0-brightgreen.svg?style=flat)](https://snakemake.github.io/)

  Snakemake is the central workflow and dependency manager of V-pipe. It determines the order in which individual tools are invoked and checks that programs do not exit unexpectedly.

  Once you have conda installed, you can in turn [use it to obtain Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda)
  (This is the recommended way to install it). Snakemake will subsequently obtain all the necessary components to V-pipe.

- **[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)**
  [![Bioconda package](https://img.shields.io/conda/dn/bioconda/fastqc.svg?label=Bioconda)](https://bioconda.github.io/recipes/fastqc/README.html)

  FastQC gives an overview of the raw sequencing data. Flowcells that have been overloaded or otherwise fail during sequencing can easily be determined with FastQC.

- **[PRINSEQ](http://prinseq.sourceforge.net/)**
  [![Bioconda package](https://img.shields.io/conda/dn/bioconda/prinseq.svg?label=Bioconda)](https://bioconda.github.io/recipes/prinseq/README.html)

  Trimming and clipping of reads is performed by PRINSEQ. It is currently the most versatile raw read processor with many customization options.

- **[Vicuna](https://www.broadinstitute.org/viral-genomics/vicuna)**
  [![Bioconda package](https://img.shields.io/conda/dn/bioconda/mvicuna.svg?label=Bioconda)](https://bioconda.github.io/recipes/mvicuna/README.html)

  Vicuna is a *de novo* assembler designed for generating rough reference contigs of viral NGS data. It can deal with the inherent heterogeneity such as high single-base heterogeneity and structural variants.

- **[InDelFixer](https://github.com/cbg-ethz/InDelFixer)**
  [![Bioconda package](https://img.shields.io/conda/dn/bioconda/indelfixer.svg?label=Bioconda)](https://bioconda.github.io/recipes/indelfixer/README.html)

  InDelFixer is a sensitive aligner employing a full Smith-Waterman alignment against a reference, used to polish up consensus.

- **[ConsensusFixer](https://github.com/cbg-ethz/ConsensusFixer)**
  [![Bioconda package](https://img.shields.io/conda/dn/bioconda/consensusfixer.svg?label=Bioconda)](https://bioconda.github.io/recipes/consensusfixer/README.html)

  ConsensusFixer is also used to polish up consensus. It computes a consensus sequence with wobbles, ambiguous bases, and in-frame insertions, from a NGS read alignment.

- **[ngshmmalign](https://github.com/cbg-ethz/ngshmmalign)**
  [![Bioconda package](https://img.shields.io/conda/dn/bioconda/ngshmmalign.svg?label=Bioconda)](https://bioconda.github.io/recipes/ngshmmalign/README.html)

  We perform the alignment of the curated NGS data using our custom ngshmmalign that takes structural variants into account. It produces multiple consensus sequences that include either majority bases or ambiguous bases.

- **[bwa](https://github.com/lh3/bwa)**
  [![Bioconda package](https://img.shields.io/conda/dn/bioconda/bwa.svg?label=Bioconda)](https://bioconda.github.io/recipes/bwa/README.html)

  In order to detect specific cross-contaminations with other probes, the Burrows-Wheeler aligner is used. It quickly yields estimates for foreign genomic material in an experiment.
  Additionally, It can be used as an alternative aligner to ngshmmalign.

- **[MAFFT](http://mafft.cbrc.jp/alignment/software/)**
  [![Bioconda package](https://img.shields.io/conda/dn/bioconda/mafft.svg?label=Bioconda)](https://bioconda.github.io/recipes/mafft/README.html)

  To standardise multiple samples to the same reference genome (say HXB2 for HIV-1), the multiple sequence aligner MAFFT is employed. The multiple sequence alignment helps in determining regions of low conservation and thus makes standardisation of alignments more robust.

- **[Samtools and bcftools](https://www.htslib.org/)**
  [![Bioconda package](https://img.shields.io/conda/dn/bioconda/samtools.svg?label=Bioconda)](https://bioconda.github.io/recipes/samtools/README.html)
  [![Bioconda package](https://img.shields.io/conda/dn/bioconda/bcftools.svg?label=Bioconda)](https://bioconda.github.io/recipes/bcftools/README.html)

  The Swiss Army knife of alignment postprocessing and diagnostics. bcftools is also used to generate consensus sequence with indels.

- **[ShoRAH](https://github.com/cbg-ethz/shorah)**
  [![Bioconda package](https://img.shields.io/conda/dn/bioconda/shorah.svg?label=Bioconda)](https://bioconda.github.io/recipes/shorah/README.html)

  The Short Reads Assembly into Haplotypes (ShoRAH) program for inferring viral haplotypes from NGS data is used to perform local haplotype reconstruction for heterogeneous viral populations by using a Gibbs sampler.

- **[LoFreq](https://csb5.github.io/lofreq/)**
  [![Bioconda package](https://img.shields.io/conda/dn/bioconda/lofreq.svg?label=Bioconda)](https://bioconda.github.io/recipes/lofreq/README.html)

  LoFreq (version 2) is SNVs and indels caller from next-generation sequencing data, and can be used as an alternative engine for SNV calling.

- **[SAVAGE](https://bitbucket.org/jbaaijens/savage)**
  [![Bioconda package](https://img.shields.io/conda/dn/bioconda/savage.svg?label=Bioconda)](https://bioconda.github.io/recipes/savage/README.html)

  SAVAGE is a tool for viral haplotype reconstruction. It can be executed in two modes: (1) using a reference sequence, or (2) assembling viral haplotypes *de novo*. We employ the latter.

- **[Haploclique](https://github.com/cbg-ethz/haploclique)**
  [![Bioconda package](https://img.shields.io/conda/dn/bioconda/haploclique.svg?label=Bioconda)](https://bioconda.github.io/recipes/haploclique/README.html)

  Viral quasispecies assembly via maximal clique finding is used as another selectable engine for global haplotype reconstruction for heterogeneous viral populations.

- **QuasiRecomb**

  QuasiRecomb performs local and global haplotype reconstruction for heterogeneous viral populations by using a hidden Markov model.

- **[SmallGenomeUtilities](https://github.com/cbg-ethz/smallgenomeutilities)**
  [![Bioconda package](https://img.shields.io/conda/dn/bioconda/smallgenomeutilities.svg?label=Bioconda)](https://bioconda.github.io/recipes/smallgenomeutilities/README.html)

  We perform genomic liftovers to standardised reference genomes using our in-house developed python library of utilities for rewriting alignments.

- **[Samtools](https://github.com/samtools/samtools)**
  [![Bioconda package](https://img.shields.io/conda/dn/bioconda/samtools.svg?label=Bioconda)](https://bioconda.github.io/recipes/samtools/README.html)

  The Swiss Army knife of alignment postprocessing and diagnostics.

- **[picard](https://broadinstitute.github.io/picard/)**
  [![Bioconda package](https://img.shields.io/conda/dn/bioconda/picard.svg?label=Bioconda)](https://bioconda.github.io/recipes/picard/README.html)

  Java tools for working with NGS data in the BAM format
