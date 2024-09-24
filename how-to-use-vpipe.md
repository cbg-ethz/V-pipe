---
layout: default
title: How to use V-pipe - V-pipe | NGS data analysis pipeline
permalink: /how-to-use-vpipe/
---

# How to use V-pipe

## Table of Contents
- [Hardware requirements](#hardware-requirements)
- [Software requirements](#software-requirements)

---

## Hardware requirements

The actual space usage depends on datasets size. Recommendation:

- RAM >= 16GiB
- Storage >= 40 GiB


- x86 CPU recommended (Apple silicon users see [OS](#OS) requirements below)
- cores >= 4

---

## Software requirements


### OS

- Linux and Mac OS X currently supported.
- Windows 10 users can install [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/about) for compatibility with Linux as bioconda [does not support Windows directly](https://bioconda.github.io/faqs.html#what-versions-are-supported).
- Mac OS X users with Apple Silicon should use the `arch -x86_64 …` command (e.g.: `env /usr/bin/arch -x86_64 /bin/zsh --login`) to run through Rosetta as bioconda support isn't stable yet.

### Software

- mamba (e.g.: from [miniforge](https://github.com/conda-forge/miniforge)) and [snakemake](https://snakemake.readthedocs.io/en/stable/) are required
  - we strongly encourage our users to deploy the pipeline [using the quick install script](https://github.com/cbg-ethz/V-pipe/blob/master/utils/README.md#quick-installer)
  - the [first tutorial, _V-pipe Installation_](https://github.com/cbg-ethz/V-pipe/blob/master/docs/tutorial_0_install.md), demonstrates how to use it
- All dependencies of V-pipe are automatically downloaded and installed by snakemake from [bioconda](https://bioconda.github.io/).
