---
jupyter:
  jupytext:
    cell_metadata_filter: -all
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---


# V-Pipe Tutorial

V-pipe is a workflow designed for the analysis of next generation sequencing (NGS) data from viral pathogens. It produces a number of results in a curated format (e.g., consensus sequences, SNV calls, local/global haplotypes). V-pipe is written using the Snakemake workflow management system.

## Requirements

V-pipe is optimized for Linux or Mac OS systems. Therefore, we recommend users with a Windows system to install WSL2 - this is not a full virtual machine but rather a way to run Windows and Linux cooperatively at the same time.  


## Organizing Data

V-Pipe takes as an input raw data in FASTQ format and depending on the user-defined configuration will output consensus sequences, SNV calls and local/global haplotypes.

V-pipe expects the input samples to be organized in a two-level hierarchy:

At the first level, input files are grouped by samples (e.g.: patients or biological replicates of an experiment).
At the second level, different datasets belonging to the same sample (e.g., from sample dates) are distinguished.
Inside the 2nd-level directory, the sub-directory `raw_data` holds the sequencing data in FASTQ format (optionally compressed with GZip).
Paired-ended reads need to be in split files with suffixes `_R1` and `_R2`.

```
samples
|───patient1
│   └───date1
│       └───raw_data
│           |───reads_R1.fastq
│           └───reads_R2.fastq
└───patient2
    |───date1
    |   └───raw_data
    |       |───reads_R1.fastq
    |       └───reads_R2.fastq
    └───date2
        └───raw_data
            |───reads_R1.fastq
            └───reads_R2.fastq
```

## Preparing a small dataset

In the directory `example_HIV_data` you find a small test dataset that you can run on your workstation or laptop.
The files will have the following structure:

```
samples
|└───CAP217
│   └───4390
│       └───raw_data
│           |───reads_R1.fastq
│           └───reads_R2.fastq
└───CAP188
    |───4
    |   └───raw_data
    |       |───reads_R1.fastq
    |       └───reads_R2.fastq
    └───30
        └───raw_data
            |───reads_R1.fastq
            └───reads_R2.fastq
```

## Install V-pipe

V-pipe uses the [Bioconda](https://bioconda.github.io/) bioinformatics software repository for all its pipeline components. The pipeline itself is implemented using [Snakemake](https://snakemake.readthedocs.io/en/stable/).

For advanced users: If your are fluent with these tools, you can:

* directly download and install [bioconda](https://bioconda.github.io/user/install.html) and [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda),
* specifiy your V-pipe configuration, and start using V-pipe

Use `--use-conda` to [automatically download and install](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management) any further pipeline dependencies. Please refer to the documentation for additional instructions.

In this present tutorial you will learn how to setup a workflow for the example dataset.

To deploy V-pipe, you can use the installation script with the following parameters:

```bash
curl -O 'https://raw.githubusercontent.com/cbg-ethz/V-pipe/master/utils/quick_install.sh'
bash quick_install.sh -p testing -w work
```

Note that

* using `-p` specifies the subdirectory where to download and install snakemake and V-pipe
* using `-w` will create a working directory and populate it. It will colloquial the references and the default `config/config.yaml`, and create a handy `vpipe` short-cut script to invoke `snakemake`.


If you get `zsh: permission denied: ./quick_install.sh`, run `chmod +x quick_install.sh` this gives the necessary permissions.

Tip: To create and populate other new working directories, you can call init_project.sh from within the new directory:

```console
mkdir -p working_2
cd working_2
../V-pipe/init_project.sh
```


## Preparation

Copy the samples directory you created in the step "Preparing a small dataset" to this working directory. You can display the directory structure with `tree sample`s or `find samples`.

```bash
mv ./samples ./testing/work/resources/
```

### Reference
If you have a reference sequences that you would like to use for read mapping and alignment, then add it to the `resources/reference/ref.fasta` directory. In our case, however, we will use the reference sequence HXB2 already provided by V-Pipe `V-pipe/resources/hiv/HXB2.fasta`.

### Preparing V-pipe's configuration

In the `work`  directory you can find the file `config.yaml`. This is where the V-Pipe configuation should be specified. See [here] (https://github.com/cbg-ethz/V-pipe/tree/master/config#readme) for the documentation of the configuration. In this tutorial we are building our own configuration therefore `virus_base_config` will remain empty. Since we are working with HIV-1, V-Pipe is providing meta information that will be used for visualisation (metainfo_file and gff_directory).

```bash
general:
    virus_base_config: ''
    aligner: "bwa"
    snv_caller: "shorah"
    haplotype_reconstruction: "haploclique"

input:
    reference: "{VPIPE_BASEDIR}/../resources/hiv/HXB2.fasta"
    metainfo_file: "{VPIPE_BASEDIR}/../resources/hiv/metainfo.yaml"
    gff_directory: "{VPIPE_BASEDIR}/../resources/hiv/gffs/"
    datadir: "{VPIPE_BASEDIR}/../../work/resources/samples"
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

Note: A YAML files use spaces as indentation, you can use 2 or 4 spaces for indentation, but no tab. There are also online YAML file validators that you might want to use if your YAML file is wrongly formatted.

## Running V-pipe


Before running check what will be executed:

```bash
cd ./testing/work/
./vpipe --dryrun
```

As this is your first run of V-pipe, it will also generate the sample collection table. Check `samples.tsv` in your editor.

Note that the samples you have downloaded have reads of length 301 only. V-pipe’s default parameters are optimized for reads of length 250. To adapt to the read length, add a third column in the tab-separated file as follows:

```bash
cat ./testing/work/samples.tsv
CAP217	4390	301
CAP188	4	301
CAP188	30	301
```

Always check the content of the `samples.tsv` file.

If you did not use the correct directory structure, this file might end up empty or some entries might be missing.
You can safely delete it and re-run with option `--dry-run` to regenerate it.

Finally, we can run the V-pipe analysis (the necessary dependencies will be downloaded and installed in conda environments managed by snakemake):

```bash
cd ./testing/work/
./vpipe -p --cores 2
```


## Output

The Wiki contains an overview of the output files. The output of the SNV calling step is aggregated in a standard [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) file, located in `samples/​{hierarchy}​/variants/SNVs/snvs.vcf`. You can open it with your favorite VCF tools for visualisation or downstream processing. It is also available in a tabular format in `samples/​{hierarchy}​/variants/SNVs/snvs.csv`.


## Swapping component

The default configuration uses ShoRAH to call the SNVs and to reconstruct the local (windowed) haplotypes.

Components of the pipeline can be swapped simply by changing the `config.yaml` file. For example to call SNVs using lofreq instead of ShoRAH use

```yaml
general:
  snv_caller: lofreq
```
