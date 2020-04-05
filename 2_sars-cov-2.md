---
layout: page
title: SARS-CoV-2
permalink: /sars-cov-2/
---

[![SARS-CoV-2]({{ "/img/sars-cov-2_header.png" | relative_url }})](https://en.wikipedia.org/wiki/Severe_acute_respiratory_syndrome_coronavirus_2)

V-pipe is part of the [SIB resources supporting SARS-CoV-2 research](https://www.sib.swiss/about-sib/news/10660-sib-resources-supporting-sars-cov-2-research).


A specifically adapted version of V-pipe to analyze high-throughput sequencing data of the SARS-CoV-2 virus has been released: 

 - Available in the [sars-cov2 branch](https://github.com/cbg-ethz/V-pipe/tree/sars-cov2) of the Github repository. 
 - Uses the NCBI reference sequence [NC_045512](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512), identical to GenBank entry [MN908947](https://www.ncbi.nlm.nih.gov/nuccore/MN908947).
 - [Default configuration](https://github.com/cbg-ethz/V-pipe/blob/sars-cov2/vpipe.config) file ready to perform SARS-CoV-2 NGS analysis.
 - New version [ShoRAH 2](https://github.com/cbg-ethz/shorah/releases) used for SNPs: with improved performance and standard [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) output.

This version will be further updated with improved reporting and results visualization features.

[Subscribe to our newsletter](https://sympa.ethz.ch/sympa/subscribe/v-pipe-users) to be informed about future updates.


## Tutorials

Check out the [tutorial - currently demonstrating the SARS-CoV-2 version]({{ "/tutorial/sars-cov2/" | relative_url }}).
