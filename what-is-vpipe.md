---
layout: default
title: What is V-pipe - V-pipe | NGS data analysis pipeline
permalink: /what-is-vpipe/
---

# What is V-pipe?

## Table of Contents
- [Introduction](#introduction)
- [What is V-pipe?](#what-is-v-pipe)
- [What can V-pipe do with my data?](#what-can-v-pipe-do-with-my-data)
- [How does V-pipe work?](#how-does-v-pipe-work)
- [Can I build my own pipeline?](#can-i-build-my-own-pipeline)
- [How can I use V-pipe?](#how-can-i-use-v-pipe)

---


## Introduction

Virus populations exist as heterogeneous ensembles of genomes within their hosts in clinical samples or among hosts in environmental samples.
This genetic diversity is associated with viral pathogenesis, virulence, and disease progression, and it can be probed using high-throughput sequencing technologies.

----

## What is V-pipe?

V-pipe is a bioinformatics pipeline that integrates various computational tools for the analysis of viral high-throughput sequencing data. 
It supports the reproducible analysis of genomic diversity in intra-host virus populations, which is often involved in viral pathogenesis and virulence.
V-pipe is targeted at researchers and professionals in patient-focused clinical applications, epidemiological surveillance, and public health management.
It is specialized in samples with mixtures of viral populations, and is positioned as a flexible, adaptable, and comprehensive tool, capable of handling different viruses and applications.

---

## What can V-pipe do with my data?

V-pipe takes as input read data obtained from a viral sequencing experiment and produces, in a single execution of the pipeline, various output files covering quality control, read alignment, and inference of viral genomic diversity on the level of both single-nucleotide variants and viral haplotypes.

---

## How does V-pipe work?

V-pipe uses the workflow management system Snakemake to determine the order in which the steps of the specified pipeline are executed and checks that the output files are produced. To simplify installation of all components, conda environments are provided. 

---

## Can I build my own pipeline?

V-pipe has a modular and extensible architecture. Users can design their own fully reproducible and transparent pipelines. Developers can test their own tools in a defined environment and contribute to the establishment of best practices for virus research and clinical diagnostics.

---

## How can I use V-pipe?

V-pipe is freely available for download:
[Get started with V-pipe]({{ site.baseurl }}/getting-started).

Further details can be found in the [documentation]({{ site.baseurl }}/documentation/) section and in the
[readme](https://github.com/cbg-ethz/V-pipe/blob/master/README.md#usage) of the software.

Should you have any further question, please do not hesitate to [contact us]({{ site.baseurl }}/contact/).
