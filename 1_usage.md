---
layout: page
title: Usage
permalink: /usage/
---

## Tutorials

- Tutorials are now directly available in V-pipe's [docs directory](https://github.com/cbg-ethz/V-pipe/blob/master/docs/README.md)
  - [tutorial_hiv.md](https://github.com/cbg-ethz/V-pipe/blob/master/docs/tutorial_hiv.md): uses HIV test data
  - [tutorial_sarscov2.md](https://github.com/cbg-ethz/V-pipe/blob/master/docs/tutorial_sarscov2.md): uses SARS-CoV-2 data from a publication

- Older tutorials built around V-pipe 2.0 and its branches are avilable here:
  - [How to use V-pipe (specifically for SARS-CoV-2 data)]({{ "/tutorial/sars-cov2/" | relative_url }}).
  - [Webinar: Applying V-pipe to SARS-Coronavirus-2 data](https://youtu.be/pIby1UooK94).

## Overview

V-pipe is designed with hierarchically organised data in mind:

```text
📁samples
├──📁patient1
│  ├──📁20100113
│  └──📁20110202
└──📁patient2
   └──📁20081130
```

Here, we have two samples from patient 1 and one sample from patient 2. All sample names should be unique such later mixups of different timepoints can be avoided.

V-pipe's parameters for the number of cores to use and the maximum memory is [specified in the config file](https://github.com/cbg-ethz/V-pipe/blob/master/config/README.md)
`config.yaml`, for instance:

```yaml
hmm_align:
  threads: 24
  leave_tmp: true
```

This instructs the `ngshmmalign` step to use 24 cores and leave the MSA temp files, which might be useful for debugging certain genomic regions.

To invoke V-pipe on the current sample set, first perform a verbose dry-run:

```bash
./vpipe -n -p
```

and after confirming that all targets are as you would expect them, perform the real run:

```bash
./vpipe -j 1
```

## Documentation

You can find more resources about using V-pipe on the [project's readme](https://github.com/cbg-ethz/V-pipe/blob/master/README.md).

### Installation

See [usage](https://github.com/cbg-ethz/V-pipe/blob/master/README.md#usage) for instructions regarding initial setup.
We strongly encourage you to deploy it [using the quick install script](https://github.com/cbg-ethz/V-pipe/blob/master/utils/README.md#quick-installer), as this is our preferred method.

### Configure options

The [config's readme](https://github.com/cbg-ethz/V-pipe/blob/master/config/README.md) gives an introduction about configuring V-pipe.
In your local installation, the file `config/config.html` contains an exhaustive list of all options that can be set in V-pipe's config file `config.yaml`.

### V-pipe as a benchmark tool

A new benchmarking is being integrated into V-pipe ahead of version 3.0. You can find a preview of this up-coming feature in the [readme of the _benchmark_ auxiliary workflow](https://github.com/cbg-ethz/V-pipe/tree/master/resources/auxiliary_workflows/benchmark/README.md). The older benchmarking provided by V-pipe 2.0 is still documented on the [old wiki](https://github.com/cbg-ethz/V-pipe/wiki/benchmark).

### Snakemake rules

V-Pipe uses [Snakemake](https://github.com/cbg-ethz/V-pipe/wiki/snakemake), a robust workflow management system, and it is possible for users to easily customize the workflow by adding or excluding rules according to their specific requirements.
