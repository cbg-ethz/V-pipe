---
layout: page
title: SARS-CoV-2
permalink: /sars-cov-2/
---

[![SARS-CoV-2]({{ "/img/sars-cov-2_header.png" | relative_url }})](https://en.wikipedia.org/wiki/Severe_acute_respiratory_syndrome_coronavirus_2)


We have released a new version of V-pipe specifically adapted  to analyze high-throughput sequencing data of the SARS-CoV-2 virus.

 - available in the [sars-cov2 branch](https://github.com/cbg-ethz/V-pipe/tree/sars-cov2) of the Github repository
 - default reference sequence is [NC_045512](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512), identical to GenBank entry [MN908947](https://www.ncbi.nlm.nih.gov/nuccore/MN908947)
 - [default configuration](https://github.com/cbg-ethz/V-pipe/blob/sars-cov2/vpipe.config) file assumes SARS-CoV-2 NGS data
 - new version [ShoRAH 2](https://github.com/cbg-ethz/shorah/releases) integrated for SNVs calling
 - all SNV calls are also reported as standard [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) files


Technical requirements:

 - Hardware: recommended RAM \>= 16GiB, Storage \>= 40 GiB
   (actual space usage depends on datasets size).
 - OS: Linux and Mac OS X currently supported
   (Windows 10 users can install [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/about)
   as bioconda [does not support Windows](https://bioconda.github.io/user/install.html) directly).
 - Software: [miniconda](https://docs.conda.io/en/latest/miniconda.html)
   and [snakemake](https://snakemake.readthedocs.io/en/stable/) required
   (see [tutorial]({{ "/tutorial/sars-cov2/#miniconda" | relative_url }})).
   All dependencies are automatically downloaded and installed by snakemake from [bioconda](https://bioconda.github.io/).


The **tutorial** demonstrates [how to use V-pipe (specifically for SARS-CoV-2 data)]({{ "/tutorial/sars-cov2/" | relative_url }}).


The V-pipe *sars-cov2* branch is continuously updated with improvements and extensions,
including on the visualization and reporting of results.

[Subscribe to our newsletter](https://sympa.ethz.ch/sympa/subscribe/v-pipe-users) to be informed about future updates.


V-pipe is part of the [SIB resources supporting SARS-CoV-2 research](https://www.sib.swiss/about-sib/news/10660-sib-resources-supporting-sars-cov-2-research).
