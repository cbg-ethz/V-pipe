---
layout: default
title: Documentation - V-pipe | NGS data analysis pipeline
nav_items:
  - name: Home
    link: "/index"
  - name: Documentation
    link: "/documentation"
    sub_items:
      - name: Tutorials
        link: "/documentation#tutorials"
      - name: General information
        link: "/documentation#general-information"
      - name: Configuration Files
        link: "/documentation#configuration-files"
      - name: How to Benchmark Tools?
        link: "/documentation#how-to-benchmark-tools"
  - name: Literature
    link: "/literature"
    sub_items:
      - name: How to Cite Us
        link: "/literature#how-to-cite-us"
      - name: Use of V-pipe
        link: "/literature#use-of-v-pipe"
      - name: V-pipe vs Competitors
        link: "/literature#v-pipe-vs-competitors"
  - name: About
    link: "/about"
    sub_items:
      - name: Funding
        link: "/about#funding"
      - name: Team
        link: "/about#team"
      - name: License and privacy policy Policy
        link: "/about#license-privacy-policy"
  - name: Contact
    link: "/contact"
    sub_items:
      - name: Chat with Us
        link: "/contact#chat-with-us"
      - name: Submit an Issue
        link: "/contact#submit-an-issue"
      - name: E-mail Us
        link: "/contact#e-mail-us"
      - name: Subscribe to our Mailing List
        link: "/contact#subscribe-to-mailing-list"
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

Due to automated texting, each copy-pastable block begins with a command entering the directory and ends with on leaving the directory:

```bash
cd tutorial/work/
# do something
cd ../..
```
Of course you don't necessarily need to do that.  You can simply remain in the working directory.

When editing files like `config.yaml`, you can use your favorite editor (`vim`, `emacs`, `nano`, [butterflies](https://xkcd.com/378/), etc.). By default our tutorials use a [_heredoc_](https://en.wikipedia.org/wiki/Here_document) to make it easier to copy-paste the blocks into bash:

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
