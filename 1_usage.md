---
layout: page
title: Usage
permalink: /usage/
---

## Tutorials

Check out the [tutorial - currently demonstrating the SARS-CoV-2 version]({{ "/tutorial/sars-cov2/" | relative_url }})

## Overview

V-pipe is designed with hierarchically organised data in mind:

```
samples
├── patient1
│   ├── 20100113
│   └── 20110202
└── patient2
    └── 20081130
```

Here, we have two samples from patient 1 and one sample from patient 2. All sample names should be unique such later mixups of different timepoints can be avoided.

V-pipe's parameters for the number of cores to use and the maximum memory is [specified in the config file](https://github.com/cbg-ethz/V-pipe/wiki/options)
`vpipe.config`, for instance:

```
[ngshmmalign]
number_cores = 24
leave_tmp = true
```

This instructs the `ngshmmalign` step to use 24 cores and leave the MSA temp files, which might be useful for debugging certain genomic regions.

To invoke V-pipe on the current sample set, first perform a verbose dry-run:

```
snakemake -n -p -s vpipe.snake
```

and after confirming that all targets are as you would expect them, perform the real run:

```
snakemake -s vpipe.snake
```

## Documentation

You can find more ressources about using V-pipe on the project's [wiki](https://github.com/cbg-ethz/V-pipe/wiki).

### Installation

See [Getting started](https://github.com/cbg-ethz/V-pipe/wiki/getting_started) for instructions regarding initial setup.

### Configure options

[V-pipe: user configurable options](https://github.com/cbg-ethz/V-pipe/wiki/options) contains a list of options that can be set in V-pipe's config file `vpipe.config`.

### V-pipe as a benchmark tool

V-pipe also provides an [unified benchmarking platform](https://github.com/cbg-ethz/V-pipe/wiki/benchmark), by incorporating two additional modules: a read simulator and a module to evaluate the accuracy of the results.

### Snakemake rules

V-Pipe uses [Snakemake](https://github.com/cbg-ethz/V-pipe/wiki/snakemake), a robust workflow management system, and it is possible for users to easily customize the workflow by adding or excluding rules according to their specific requirements.
