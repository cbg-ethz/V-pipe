---
layout: page
title: Pipeline overview
permalink: /pipeline/
---


## Pipeline overview

![Pipeline](../img/pipeline.svg)


## Dependencies

- **[Snakemake](https://snakemake.bitbucket.io)**
  [![Bioconda package](https://img.shields.io/conda/dn/bioconda/snakemake.svg?label=Bioconda)](http://bioconda.github.io/recipes/snakemake/README.html)
  [![Snakemake](https://img.shields.io/badge/snakemake-≥4.8.0-brightgreen.svg?style=flat)](https://snakemake.bitbucket.io)

  Snakemake is the central workflow and dependency manager of V-pipe. It determines the order in which individual tools are invoked and checks that programs do not exit unexpectedly.

- **[Conda](https://conda.io/docs/index.html)**

  Conda is an open source package management system and environment management system. V-pipe uses it to automatically obtain reproducible environments and simplify installation of the individual components of the pipeline, thanks to the [Bioconda channel](https://bioconda.github.io) - a distribution of bioinformatics software.

  See [Documentation](http://conda.io/docs/install/quick.html) to install it.

- **FastQC**

  FastQC gives an overview of the raw sequencing data. Flowcells that have been overloaded or otherwise fail during sequencing can easily be determined with FastQC.

- **[PRINSEQ](http://prinseq.sourceforge.net/)**
  [![Bioconda package](https://img.shields.io/conda/dn/bioconda/prinseq.svg?label=Bioconda)](http://bioconda.github.io/recipes/prinseq/README.html)

  Trimming and clipping of reads is performed by PRINSEQ. It is currently the most versatile raw read processor with many customization options.

- **[Vicuna](https://www.broadinstitute.org/viral-genomics/vicuna)**
  [![Bioconda package](https://img.shields.io/conda/dn/bioconda/mvicuna.svg?label=Bioconda)](https://bioconda.github.io/recipes/mvicuna/README.html)

  Vicuna is a de novo assembler designed for generating rough reference contigs of viral NGS data. It can deal with the inherent heterogeneity such as high single-base heterogeneity and structural variants.

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
  It can also optionnally be used as an alternative to ngshmmalign.

- **[MAFFT](http://mafft.cbrc.jp/alignment/software/)**
  [![Bioconda package](https://img.shields.io/conda/dn/bioconda/mafft.svg?label=Bioconda)](https://bioconda.github.io/recipes/mafft/README.html)

  To standardise multiple samples to the same reference genome (say HXB2 for HIV-1), the multiple sequence aligner MAFFT is employed. The multiple sequence alignment helps in determining regions of low conservation and thus makes standardisation of alignments more robust.

- **[ShoRAH](https://github.com/cbg-ethz/shorah)**
  [![Bioconda package](https://img.shields.io/conda/dn/bioconda/shorah.svg?label=Bioconda)](https://bioconda.github.io/recipes/shorah/README.html)

  The Short Reads Assembly into Haplotypes (ShoRAH) program for inferring viral haplotypes from NGS data is used to perform local haplotype reconstruction for heterogeneous viral populations by using a Gibbs sampler.

- **[LoFreq](https://csb5.github.io/lofreq/)**
  [![Bioconda package](https://img.shields.io/conda/dn/bioconda/lofreq.svg?label=Bioconda)](https://bioconda.github.io/recipes/lofreq/README.html)

  LoFreq version 2 is a fast and sensitive variant-caller for inferring SNVs and indels from next-generation sequencing data, and can be optionnally used as an alternative engine for SNP calling.

- **[Savage](https://bitbucket.org/jbaaijens/savage)**
  [![Bioconda package](https://img.shields.io/conda/dn/bioconda/savage.svg?label=Bioconda)](https://bioconda.github.io/recipes/savage/README.html)

  SAVAGE is a computational tool for reconstructing individual haplotypes of intra-host virus strains (a viral quasispecies) without the need for a high quality reference genome. SAVAGE makes use of either FM-index based data structures or ad-hoc consensus reference sequence for constructing overlap graphs from patient sample data. It is one of the selectable engines for global haplotype reconstruction for heterogeneous viral populations.

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
