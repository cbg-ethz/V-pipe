![Logo](https://cbg-ethz.github.io/V-pipe/img/logo.svg)

V-pipe is a workflow designed for analysis of next generation sequencing (NGS) data from viral pathogens. It produces a number of results in a curated format.

[![Snakemake](https://img.shields.io/badge/snakemake-≥4.8.0-blue.svg?style=flat-square)](https://snakemake.bitbucket.io)
[![bio.tools](https://img.shields.io/badge/bio-tools-blue.svg?style=flat-square)](https://bio.tools/V-Pipe)
[![License: Apache-2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg?style=flat-square)](https://opensource.org/licenses/Apache-2.0)

## Quick start

Instructions to type in a shell

1. [Install](https://docs.conda.io/projects/continuumio-conda/en/latest/user-guide/install/index.html) miniconda3

### Linux

  To obtain the installer for linux use the following:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

  Then, install miniconda,

```
sh Miniconda3-latest-Linux-x86_64.sh
```

### MacOS

  To obtain the installer for MacOS, you can [download](https://docs.conda.io/en/latest/miniconda.html) it manually or use wget:
```
curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
```

  Then, install miniconda,

```
sh Miniconda3-latest-MacOSX-x86_64.sh
```

2. Create conda virtual environment

```
conda create -n V-pipe -c bioconda python=3.8 snakemake-minimal=5.14.0
conda activate V-pipe
```

Make sure to use `source activate V-pipe` everytime you want to run V-pipe

3. Get V-pipe

```
git clone https://github.com/cbg-ethz/V-pipe.git /path/to/V-pipe
```

## Running V-pipe

First, open a terminal and change into the **working directory** where input files are stored (i.e., the reference and the sequencing reads). We use a [two-level](https://github.com/cbg-ethz/V-pipe/wiki/getting-started#input-files) directory hierarchy and we expect sequencing reads in a folder name `raw_data`. To initialize a project,

```
/path/to/V-pipe/init_project.sh
```

Before actually running the pipeline, we advise to check whether output files can be created from the inputs, using the `--dryrun` option.
```
./vpipe --dryrun
```

Further details can be found in the [wiki](https://github.com/cbg-ethz/V-pipe/wiki) pages.

## Dependencies

- **conda**

  Conda is a cross-platform package management system and an environment manager application.

- **Snakemake**

  Snakemake is the central workflow and dependency manager of V-pipe. It determines the order in which individual tools are invoked and checks that programs do not exit unexpectedly.

- **VICUNA**

  VICUNA is a *de novo* assembly software designed for populations with high mutation rates. It is used to build an initial reference for mapping reads with ngshmmalign aligner when a `references/cohort_consensus.fasta` file is not provided. Further details can be found in the [wiki](https://github.com/cbg-ethz/V-pipe/wiki/getting-started#input-files) pages.

### Computational tools 
Other dependencies are managed by using isolated conda environments per rule, and below we list some of the computational tools integrated in V-pipe:

- **PRINSEQ**

  Trimming and clipping of reads is performed by PRINSEQ. It is currently the most versatile raw read processor with many customization options.

- **Vicuna**

  Vicuna is a de novo assembler designed for generating rough reference contigs of viral NGS data. It can deal with the inherent heterogeneity such as high single-base heterogeneity and structural variants.

- **ngshmmalign**

  We perform the alignment of the curated NGS data using our custom ngshmmalign that takes structural variants into account. It produces multiple consensus sequences that include either majority bases or ambiguous bases.

- **bwa**

  In order to detect specific cross-contaminations with other probes, the Burrows-Wheeler aligner is used. It quickly yields estimates for foreign genomic material in an experiment.

- **MAFFT**

  To standardise multiple samples to the same reference genome (say HXB2 for HIV-1), the multiple sequence aligner MAFFT is employed. The multiple sequence alignment helps in determining regions of low conservation and thus makes standardisation of alignments more robust.

- **Samtools**

  The Swiss Army knife of alignment postprocessing and diagnostics.

- **SmallGenomeUtilities**

  We perform genomic liftovers to standardised reference genomes using our in-house developed python library of utilities for rewriting alignments.

- **ShoRAH**

  ShoRAh performs SNV calling and local haplotype reconstruction by using bayesian clustering.

- **HaploClique and SAVAGE**

  We use HaploClique or SAVAGE to perform global haplotype reconstruction for heterogeneous viral populations by using an overlap graph.

## Contributions

- [Susana Posada Céspedes](https://orcid.org/0000-0002-7459-8186)
- [David Seifert](https://orcid.org/0000-0003-4739-5110)
- Tobias Marschall
- [Niko Beerenwinkel](https://orcid.org/0000-0002-0573-6119)

## Contact

We encourage users to use the [issue tracker](https://github.com/cbg-ethz/V-pipe/issues). For further enquiries, you can also contact the V-pipe Dev Team <v-pipe@bsse.ethz.ch>.
