---
layout: page
title: Usage
permalink: /usage/
---

## Tutorials

- [How to use V-pipe (specifically for SARS-CoV-2 data)]({{ "/tutorial/sars-cov2/" | relative_url }}).
- [Webinar: Applying V-pipe to SARS-Coronavirus-2 data](https://youtu.be/pIby1UooK94).

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

V-pipe's parameters for the number of cores to use and the maximum memory is [specified in the config file](https://github.com/cbg-ethz/V-pipe/blob/master/config/README.md)
`config.yaml`, for instance:

```yaml
ngshmmalign:
  number_cores: 24
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

V-pipe also provides an [unified benchmarking platform](https://github.com/cbg-ethz/V-pipe/wiki/benchmark), by incorporating two additional modules: a read simulator and a module to evaluate the accuracy of the results.

### Snakemake rules

V-Pipe uses [Snakemake](https://github.com/cbg-ethz/V-pipe/wiki/snakemake), a robust workflow management system, and it is possible for users to easily customize the workflow by adding or excluding rules according to their specific requirements.
