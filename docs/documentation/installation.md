---
jupyter:
  jupytext:
    cell_metadata_filter: -all
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.5
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

<!-- markdownlint-configure-file { "MD010": { "ignore_code_languages" : [ "tsv", "bash" ] } } -->
# Installation

The present tutorial will show you how to install V-pipe and the dependencies required to start using it - bioconda, conda-froge mamba and snakemake - before continuing with other tutorials and analysing virus data.

## For the impatient

Download the install script and run it with the following parameters:

```bash
curl -O 'https://raw.githubusercontent.com/cbg-ethz/V-pipe/master/utils/quick_install.sh'
bash quick_install.sh -p vp-analysis -w work
```

## Requirements

V-pipe is optimized for Linux or Mac OS systems, and we heavily rely on bioconda, which isn't supported on Windows. Therefore, we recommend users with a Windows system to [install WSL2](https://learn.microsoft.com/en-us/windows/wsl/install).

## Quick install V-pipe and conda

V-pipe uses the [Bioconda](https://bioconda.github.io/) bioinformatics software repository for all its pipeline components. The pipeline itself is implemented using [Snakemake](https://snakemake.readthedocs.io/en/stable/). Although you can install all the dependencies manually, we recommend using our install quick install script:

```bash
curl -O 'https://raw.githubusercontent.com/cbg-ethz/V-pipe/master/utils/quick_install.sh'
bash quick_install.sh -p vp-analysis -w work
```

The script `quick_install.sh` has the following options:

* using `-p` specifies the subdirectory where to download and install snakemake and V-pipe
* using `-w` will create a working directory and populate it. It will create a boilerplate `config/config.yaml`, and create a handy `vpipe` short-cut script to invoke `snakemake`.
* an additional option `-b` (not demonstrated above) allows to install a spefic branch or tagged version. If nothing is specified, the master branch will be installed.

```{tip}
To create and populate other new working directories, you can call `init_project.sh` from within the new directory:

```bash
cd vp-analysis/

mkdir -p working_2
cd working_2
../V-pipe/init_project.sh

```

After running the `quick_install.sh` script, you should have a directory structure like this:

```text
vp-analysis
├── Mambaforge-Darwin-x86_64.sh
├── V-pipe # cloned from https://github.com/cbg-ethz/V-pipe
│   ├── CONTRIBUTING.md
│   └── ..
├── mambaforge # installation of dependencies including snakemake
│   ├── LICENSE.txt
│   └── ..
└── work # working directory
    ├── config.yaml
    └── vpipe
```

- `vp-analysis` is the main directory where we will store everything.
- `mambaforge` is the directory where conda will be installed including the dependencies to start using V-pipe.
- `V-pipe` is the directory where V-pipe's code will be downloaded from GitHub
- `work` finally, each analysis of virus data will be performed in a directory like `work…`. If you start a new analysis of a dataset, you can create a new directory, run `init_project.sh` inside the directory and get started.

## Other installation options

### Cloning the repository 

The V-pipe repository contains a snakemake pipeline. In order to run it directly with snakemake, clone the repository with:

```sh
git clone https://github.com/cbg-ethz/V-pipe.git
```

If you haven't already done so, install snakemake by using the [official instructions](https://github.com/cbg-ethz/V-pipe.git), and you can run the pipeline with `snakemake --use-conda`. 

### Using Docker

```{note}
Note: the [docker image](https://github.com/cbg-ethz/V-pipe/pkgs/container/v-pipe) is only setup with components to run the workflow for HIV and SARS-CoV-2 virus base configurations.
Using V-pipe with other viruses or configurations might require internet connectivity for additional software components.
```

Create `config.yaml` and then populate the directory containing raw reads, typically `samples/`.
For example, the following config file could be used:

```yaml
general:
  virus_base_config: hiv

output:
  snv: true
  local: true
  global: false
  visualization: true
  QA: true
```

Then execute:

```bash
docker run --rm -it -v $PWD:/work ghcr.io/cbg-ethz/v-pipe:master --jobs 4 --printshellcmds --dry-run
```

### Using Snakedeploy

Install snakedeploy according to the [official instructions](https://snakedeploy.readthedocs.io/en/latest/getting_started/installation.html).

Snakemake's [official workflow installer Snakedeploy](https://snakemake.github.io/snakemake-workflow-catalog/?usage=cbg-ethz/V-pipe) can now be used:

```bash
snakedeploy deploy-workflow https://github.com/cbg-ethz/V-pipe --tag master .
# edit config/config.yaml and provide samples/ directory
snakemake --use-conda --jobs 4 --printshellcmds --dry-run
```
