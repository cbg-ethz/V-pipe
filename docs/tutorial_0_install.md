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
# V-Pipe Installation

V-pipe is a workflow designed for the analysis of next-generation sequencing (NGS) data from viral pathogens. It produces a number of results in a curated format (e.g., consensus sequences, SNV calls, local/global haplotypes). V-pipe is written using the Snakemake workflow management system.

The present tutorial will show you how to install V-pipe and the dependencies required to start using it - bioconda, conda-froge mamba and snakemake - before continuing with other tutorials and analysing virus data.

## Requirements

V-pipe is optimized for Linux or Mac OS systems, and bioconda isn't supported on Windows. Therefore, we recommend users with a Windows system to [install WSL2](https://learn.microsoft.com/en-us/windows/wsl/install) - this is not a full virtual machine but rather a way to run Windows and Linux cooperatively at the same time.


## Organizing Software

We will organise our software in the following tree structure, which will be reused in all subsequent tutorials:

```text
📁 [HOME]
└───📁vp-analysis
    ├───📁V-pipe      # V-pipe checked out from Github
    ├───📁Miniforge3  # bioconda + conda-forge + mamba + Snakemake
    ├───📁work        # work directories
    ├───📁work-tests  #  …
    └───📁 …          #  …
```

- `vp-analysis` is the main directory where we will store everything.
- `Miniforge3` is the directory where conda will be installed including the dependencies to start using V-pipe.
- `V-pipe` is the directory where V-pipe's code will be downloaded from GitHub
- finally, each analysis of virus data will be performed in a directory like `work…`, which holds the configuration and the sequencing data for that particular analysis.


## Install V-pipe and conda from scratch

V-pipe uses the [Bioconda](https://bioconda.github.io/) bioinformatics software repository for all its pipeline components. The pipeline itself is implemented using [Snakemake](https://snakemake.readthedocs.io/en/stable/).

For advanced users: If your are fluent with these tools, see [below](#fluent-users)

In this short tutorial, you will learn how to setup a workflow for the various examples in the analysis tutorials.

To deploy V-pipe, you can use the installation script with the following parameters:

```bash
curl -O 'https://raw.githubusercontent.com/cbg-ethz/V-pipe/master/utils/quick_install.sh'
bash quick_install.sh -p vp-analysis -w work
```

> **Note**:
> * using `-p` specifies the subdirectory where to download and install snakemake and V-pipe
> * using `-w` will create a working directory and populate it. It will colloquial the references and the default `config/config.yaml`, and create a handy `vpipe` short-cut script to invoke `snakemake`.
> * an additional option `-b` (not demonstrated above) allows to install a spefic branch or tagged version. If nothing is specified, the master branch will be installed.

If you get `zsh: permission denied: ./quick_install.sh`, run `chmod +x quick_install.sh` this gives the necessary permissions.


**Tip:** To create and populate other new working directories, you can call `init_project.sh` from within the new directory:

```bash
cd vp-analysis/

mkdir -p working_2
cd working_2
../V-pipe/init_project.sh

# now edit config.yaml, samples.tsv, run your analysis etc.

cd -
```

### Analyse data

Now that you have setup the software necessary to start using V-pipe, you can follow with one of the tutorials showing you the analysis of viral sequencing data:

- [tutorial_hiv.md](tutorial_hiv.md): uses HIV test data
- [tutorial_sarscov2.md](tutorial_sarscov2.md): uses SARS-CoV-2 data from a publication


## Fluent users

For advanced users: If your are fluent with these tools, you can:

* directly download and install [Miniforge3](https://github.com/conda-forge/miniforge#Download), setup [bioconda](https://bioconda.github.io/index.html#usage) and install [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda),
* specifiy your V-pipe configuration, and start using V-pipe

Use `--use-conda` to [automatically download and install](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management) any further pipeline dependencies. Please refer to the documentation for additional instructions.


### Reusing an existing conda installation

[TO BE CONTINUED]
