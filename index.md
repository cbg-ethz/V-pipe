---
layout: default
image: img/workflow_cartoon2.svg
image: img/workflow_cartoon2.png
---


V-pipe: A bioinformatics pipeline for viral sequencing data

{% include flashNews.html urgent="0" msg="New version v2.99.3 of V-pipe <a href='https://github.com/cbg-ethz/V-pipe/releases'>has been released</a><br />includes primer trimming on cohorts with multiple different protocols" %}


# Introduction

Virus populations exist as heterogeneous ensembles of genomes within their hosts.
This genetic diversity is associated with viral pathogenesis, virulence, and disease progression, and it can be probed using high-throughput sequencing technologies.

![Cartoon](img/workflow_cartoon2.svg)

<div class="calltoaction">
  <div><a class="hrefbut" href="https://github.com/cbg-ethz/V-pipe#using-quick-install-script" style="color:black;"><img src="img/mark-github.svg" alt="GitHub" /> Install from GitHub!</a></div>
  <div><a class="hrefbut" href="https://github.com/cbg-ethz/V-pipe/pkgs/container/v-pipe" style="color:black;"><img src="img/fa-brand_fa-docker.svg" alt="Docker" /> Get the Docker image!</a></div>
  <div><a class="hrefbut" href="https://snakemake.github.io/snakemake-workflow-catalog/?usage=cbg-ethz/V-pipe" style="color:black;"><img src="img/logo-snake.svg" alt="Snakemake" /> Snakedeploy the workflow!</a></div>
</div>

[![bio.tools](https://img.shields.io/badge/bio-tools-orange.svg?style=flat)](https://bio.tools/V-Pipe)
[![ExPASy](https://img.shields.io/badge/expasy-resource-red.svg?style=flat)](https://www.expasy.org/resources/v-pipe)
[![SIB](https://img.shields.io/badge/sib-resource-red.svg?style=flat)](https://www.sib.swiss/research-infrastructure/database-software-tools/sib-resources#v-pipe)
[![workflowhub](https://img.shields.io/badge/workflow-hub-teal.svg?style=flat)](https://workflowhub.eu/workflows/301)
[![License: Apache-2.0](https://img.shields.io/badge/License-Apache_2.0-yellow.svg?style=flat)](https://opensource.org/licenses/Apache-2.0)

----

## What is V-pipe?

V-pipe is a bioinformatics pipeline that integrates various computational tools for the analysis of viral high-throughput sequencing data. It supports the reproducible analysis of genomic diversity in intra-host virus populations, which is often involved in viral pathogenesis and virulence.


## What can V-pipe do with my data?

V-pipe takes as input read data obtained from a viral sequencing experiment and produces, in a single execution of the pipeline, various output files covering quality control, read alignment, and inference of viral genomic diversity on the level of both single-nucleotide variants and viral haplotypes.


## How does V-pipe work?

V-pipe uses the workflow management system Snakemake to determine the order in which the steps of the specified pipeline are executed and checks that the output files are produced. To simplify installation of all components, conda environments are provided. 


## Can I build my own pipeline?

V-pipe has a modular and extensible architecture. Users can design their own fully reproducible and transparent pipelines. Developers can test their own tools in a defined environment and contribute to the establishment of best practices for virus research and clinical diagnostics.


## How can I use V-pipe?

V-pipe is freely available for download from [GitHub](https://github.com/cbg-ethz/V-pipe).
Further details on how to run the pipeline, as well as [test datasets](https://github.com/cbg-ethz/V-pipe/blob/master/tests/data/),
can be found in the [readme](https://github.com/cbg-ethz/V-pipe/blob/master/README.md#usage),
accessible through the V-pipe website's [Usage tab](usage/).

Should you have any further question, please do not hesitate to [contact us](contact/).

## Video presentation of V-pipe

{% include youtubePlayer.html id="qHEUVJZsgE4" %}
