---
layout: default
title: Documentation - V-pipe | NGS data analysis pipeline
permalink: /documentation/
---

# Documentation

## Table of Contents
- [Tutorials](#tutorials)
- [General Information](#general-information)
- [Configuration Files](#configuration-files)
- [How to Benchmark Tools](#how-to-benchmark-tools)

---

## Tutorials

We strongly advise our users to start discovering V-pipe by looking at the tutorials

You can find several tutorials in the `docs/` directory:

### Getting V-pipe installed

- [V-pipe Installation](https://github.com/cbg-ethz/V-pipe/blob/master/docs/tutorial_0_install.md)

### Viruses

- [V-Pipe HIV Tutorial](https://github.com/cbg-ethz/V-pipe/blob/master/docs/tutorial_hiv.md): uses HIV test data
- [SARS-CoV-2 Tutorial](https://github.com/cbg-ethz/V-pipe/blob/master/docs/tutorial_sarscov2.md): uses SARS-CoV-2 data from a publication

### Note about the tutorials

Due to automated testing, each copy-pastable block begins with a command entering the directory and ends with one leaving the directory:

```bash
cd tutorial/work/
# do something
cd ../..
```
Of course, you don't necessarily need to do that.  You can simply remain in the working directory.

When editing files like `config.yaml`, you can use your favorite editor (`vim`, `emacs`, `nano`, [butterflies](https://xkcd.com/378/), etc.).
By default, our tutorials use a [_heredoc_](https://en.wikipedia.org/wiki/Here_document) to make it easier to copy-paste the blocks into bash:

```bash
cat > config.yaml <<EOF
general:
    virus_base_config: 'sars-cov-2'
EOF
```

---

## General Information

You'll find a short introduction to V-pipe in the [README file at the root of V-pipe](https://github.com/cbg-ethz/V-pipe/blob/master/README.md#usage)

## Configuration files

In order to start using V-pipe, you need to provide three things:

1. Samples in a specific directory structure
2. (optional) TSV file listing the samples
3. Configuration file

To configure V-pipe refer to the documentation present in [config/README.md](https://github.com/cbg-ethz/V-pipe/blob/master/config/README.md).

### configuration manual

More information about all the available configuration options and an exhaustive list can be found in the file `config/config.html` in your local installation of V-pipe (temporarily not available online).


## How to benchmark

How to perform benchmarks is documented in [resources/auxiliary_workflows/benchmark/README.md](https://github.com/cbg-ethz/V-pipe/tree/master/resources/auxiliary_workflows/benchmark/README.md)

## Video tutorial

Comming soon
