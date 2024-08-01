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
# V-Pipe HIV Tutorial

The present tutorial will show you how to apply V-pipe on HIV sequencing data. The data originate from the publication Abrahams et al. (2019), Science translational medicine 11.513 (DOI: 10.1126/scitranslmed.aaw5589). For testing the installation the same data is used. In this tutorial we will show you how to set up the configuration file, run the pipeline and interpret the output.

## Requirements

The tutorial assumes that you have installed V-pipe using the [quick install installation documentation](quick-install-v-pipe-and-conda), and that the workflow is setup with the following structure:

```text
vp-analysis
├── V-pipe
├── mambaforge
└── work
```

- `vp-analysis` is the main directory where you have installed V-pipe
- `V-pipe` is the directory with V-pipe's own code
- `mambaforge` has dependencies to start using V-pipe (bioconda, conda-forge, mamba, snakemake)
- `work` is the directory where you have performed the test analysis

## Setting up the work directory

We will create a fresh work directory for this tutorial. A V-pipe work directory typically contains the following files and directories:

- `config.yaml`: the configuration file. For example to tell V-pipe where to find the samples and the reference genome. All configuration options are described in the [configuration schema](configuring-the-workflow)
- `samples.tsv`: a tab-separated file listing the samples to be processed. The first two columns are mandatory and represent the hierarchical levels of the samples. The third and fourth column are optional and contain the read length and protocol name. 
- `vpipe`: a wrapper script to start the workflow
- `samples/`: the directory containing the raw data of the samples

And after running the workflow:

- `results/`: the results of the workflow
- `.snakemake`: the directory containing the snakemake working files

For your convenience, you can set up a boilerplate working directory with the script `init.sh`. This will copy a `config.yaml` and the `vpipe` wrapper script to get started: 

```bash
cd vp-analysis

mkdir -p work_hiv
cd work_hiv
../V-pipe/init_project.sh
```

## Preparing the dataset

As described in [configuration](organizing-data) V-pipe expects the input samples to be organized in a two-level hierarchy. In the directory `vp-analysis/V-pipe/documentation/example_HIV_data` you can find a small dataset in the correct format that we will use in this tutorial. The files will have the following structure:

```text
📁samples
├───📁CAP217
│   └───📁4390
│       └───📁raw_data
│           ├───🧬reads_R1.fastq
│           └───🧬reads_R2.fastq
└───📁CAP188
    │───📁4
    │   └───📁raw_data
    │       ├───🧬reads_R1.fastq
    │       └───🧬reads_R2.fastq
    └───📁30
        └───📁raw_data
            ├───🧬reads_R1.fastq
            └───🧬reads_R2.fastq
```

Now, copy this dataset to the `work_hiv` directory:

```bash
cp -r ../V-pipe/docs/example_HIV_data/samples .
```

## Configuration

### References

V-pipe comes with pre-configured references. You can specify the species you are using in the configuration file at `general.virus_base_config`. If you are working with a reference that is not pre-configured, you can specify it in the configuration file at `input.reference`. For more information see the [documentation](configuring-the-workflow).

In this tutorial we will use a reference already available in `V-pipe/resources/hiv/HXB2.fasta`, but instead of specifing `virus_base_config` we will specify the reference directly in the configuration file. We will also seperately specify the genome annotation (gff file) and the metainfo file. 

### Populating `config.yaml`

In the `work_hiv`  directory you can find the file `config.yaml`. Open it in your editor and add the following content:

```yaml
general:
    virus_base_config: ""
    aligner: bwa
    snv_caller: shorah
    haplotype_reconstruction: haploclique

input:
    # the references are part of the repository in this case:
    reference: "../V-pipe/resources/hiv/HXB2.fasta"
    metainfo_file: "../V-pipe/resources/hiv/metainfo.yaml"
    gff_directory: "../V-pipe/resources/hiv/gffs/"
    datadir: samples/
    # we specify the read length here, as it is not the default 250:
    read_length: 301
    samples_file: samples.tsv
    paired: true

snv:
    consensus: false

output:
    snv: true
    local: true
    global: true
    visualization: true
    QA: false
    diversity: true
```

```{note}
A YAML files use spaces as indentation, you can use 2 or 4 spaces for indentation, but **no tab**. There are also [online YAML file validators](https://www.yamllint.com/) that you might want to use if your YAML file is wrongly formatted.
```

## Running V-pipe

Before running check what will be executed:

```bash
cd vp-analysis/work_hiv/

./vpipe --dryrun

```

As this is your first run of V-pipe, it will automatically generate the sample collection table (`samples.tsv`). Check `samples.tsv` in your editor. It is always a good idea check the content of the `samples.tsv` file, as it is used to collect the samples for the analysis. Of course, you can also provide `samples.tsv` yourself, before running the pipeline. If you did not use the expected directory structure, this file might end up empty or some entries might be missing. If so, you can safely delete it and re-run with option `--dry-run` to regenerate it. More information on the `samples.tsv` file can be found in the [documentation](setting-up-samplestsv).

Finally, we can run the V-pipe analysis. The first run will take a while because it will install all necessary software dependencies with conda:

```bash
cd vp-analysis/work_hiv/

./vpipe -p --cores 2
# -p and --cores (and all other options) are passed to snakemake. -p is for printing shell cmds. 
# takes a while to run, needs to install packages
```

```{note}
Note that `vpipe` is a wrapper for `snakemake`. All options that are passed to `vpipe` are options to snakemake. More information about snakemake options can be found in the [snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html).
```

## Output

The output of the SNV calling step is aggregated in a standard [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) file, located in `results/​{hierarchy}​/variants/SNVs/snvs.vcf`. You can open it with your favorite VCF tools for visualisation or downstream processing. It is also available in a tabular format in `results/​{hierarchy}​/variants/SNVs/snvs.csv`.
