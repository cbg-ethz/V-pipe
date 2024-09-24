---
layout: page
title: SARS-CoV-2
permalink: /sars-cov-2/
---

[![SARS-CoV-2]({{ "/assets/img/sars-cov-2_header.png" | relative_url }})](https://en.wikipedia.org/wiki/Severe_acute_respiratory_syndrome_coronavirus_2)


Since version v2.99.0, V-pipe integrates specific adaptations to analyze high-throughput sequencing data of SARS-CoV-2:

 - available in the [master branch](https://github.com/cbg-ethz/V-pipe/tree/master) of the Github repository and in the [releases](https://github.com/cbg-ethz/V-pipe/releases)
 - the [virus base configuration](https://github.com/cbg-ethz/V-pipe/blob/master/config/README.md#virus-base-config)
     `sars-cov-2` provides
     [defaults](https://github.com/cbg-ethz/V-pipe/blob/master/config/sars-cov-2.yaml)
     to analyze SARS-CoV-2 NGS data.
   - default reference sequence is [NC_045512](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512) (identical to GeneBank entry [MN908947](https://www.ncbi.nlm.nih.gov/nuccore/MN908947))
   - [ShoRAH 2](https://github.com/cbg-ethz/shorah/releases) integrated for SNV calling
 - all SNV calls are also reported as standard [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) files


Technical requirements:

 - Hardware: recommended RAM \>= 16GiB, Storage \>= 40 GiB
   (actual space usage depends on datasets size).
 - OS: Linux and Mac OS X currently supported
   (Windows 10 users can install [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/about)
   as bioconda [does not support Windows](https://bioconda.github.io/user/install.html) directly).
 - Software: [miniconda](https://docs.conda.io/en/latest/miniconda.html)
   and [snakemake](https://snakemake.readthedocs.io/en/stable/) required;
   we strongly encourage our users to deploy the pipeline [using the quick install script](https://github.com/cbg-ethz/V-pipe/blob/master/utils/README.md#quick-installer).
   All dependencies are automatically downloaded and installed by snakemake from [bioconda](https://bioconda.github.io/).


The **tutorial** demonstrates [how to use V-pipe (specifically for SARS-CoV-2 data)](https://github.com/cbg-ethz/V-pipe/blob/master/docs/tutorial_sarscov2.md).

V-pipe is continuously updated with improvements and extensions,
including on the visualization and reporting of results.

[Subscribe to our newsletter](https://sympa.ethz.ch/sympa/subscribe/v-pipe-users) to be informed about future updates.

V-pipe is part of the [SIB resources supporting SARS-CoV-2 research](https://www.sib.swiss/about-sib/news/10660-sib-resources-supporting-sars-cov-2-research),
used to process NGS data for the [Swiss SARS-CoV-2 Sequencing Consortium (S3C)](https://bsse.ethz.ch/cevo/research/sars-cov-2/swiss-sars-cov-2-sequencing-consortium.html)
and for the [Surveillance of SARS-CoV-2 genomic variants in wastewater](https://bsse.ethz.ch/cbg/research/computational-virology/sarscov2-variants-wastewater-surveillance.html).



*[SNV]: Single Nucleotide Variant
*[NGS]: Next Generation Sequencing
*[S3C]: Swiss SARS-CoV-2 Sequencing Consortium
