---
layout: default
---


V-pipe is a bioinformatics pipeline for the analysis of next-generation sequencing data derived from intra-host viral populations.


<div align="center" style="margin: 2%;">
  <a class="hrefbut" href="https://github.com/cbg-ethz/V-pipe" style="color:black;"><img src="img/mark-github.svg" alt="GitHub" /> Link to repository!</a>
</div>

[![bio.tools](https://img.shields.io/badge/bio-tools-orange.svg?style=flat)](https://bio.tools/V-Pipe)
[![ExPASy](https://img.shields.io/badge/expasy-resource-red.svg?style=flat)](https://www.expasy.org/resources/search/querytext:v-pipe)
[![License: Apache-2.0](https://img.shields.io/badge/License-Apache_2.0-yellow.svg?style=flat)](https://opensource.org/licenses/Apache-2.0)

----

## Introduction

Virus populations exist as heterogeneous ensembles of genomes within their hosts. 
This genetic diversity is associated with viral pathogenesis, virulence, and disease progression, and it can be probed using high-throughput sequencing technologies. 
V-pipe is a bioinformatics pipeline for the analysis of viral genomic data that allows for assessing intra-host viral genetic diversity. 
Integrating various open-source software packages, V-pipe supports quality control, read mapping, error correction, and viral haplotype reconstruction. 
It provides standardized, transparent, and reproducible workflows for research and diagnostic applications.

![Cartoon](img/workflow_cartoon.svg)

## Features
- Reference guided genome assembly, useful for phylogenetic inference for instance
- Detailed report, including quality overview, fraction of failed reads
- Contamination checking, in order to detect flowcell cross contamination from other sequencing runs
- Genomic information on three resolution scales:
  * SNV: Frequencies and positions of single nucleotide variants that differ from a control population
  * Local: co-occurrence of SNVs in regions that are as long as the average read
  * Global: haplotypes of larger segments of viral genomes
