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

V-pipe is a workflow designed for the analysis of next generation sequencing (NGS) data from viral pathogens. It produces a number of results in a curated format (e.g., consensus sequences, SNV calls, local/global haplotypes). V-pipe is written using the Snakemake workflow management system.

The present tutorial will show you how to apply V-pipe on HIV sequencing data.

## Requirements

The tutorial assumes that you have [installed V-pipe using the installation tutorial](tutorial_0_install.md), and that the workflow is setup with the following structure:

```text
📁 [HOME]
└───📁vp-analysis
    ├───📁V-pipe      # V-pipe checked out from Github
    ├───📁Miniforge3  # bioconda + conda-forge + mamba + Snakemake
    ├───📁work        # work directories
    ├───📁work-tests  #  …
    └───📁 …          #  …
```

- `vp-analysis` is the main directory where we installed everything in the previous tutorial
- `Miniforge3` has dependencies to start using V-pipe (bioconda, conda-forge, mamba, snakemake)
- `V-pipe` is the directory with V-pipe's own code
- and for this tutorial we will create a directory like `work…`, which will hold the configuration and the sequencing data for our analysis.



## Organizing Data

V-Pipe takes as an input raw data in FASTQ format and depending on the user-defined configuration will output consensus sequences, SNV calls and local/global haplotypes.

V-pipe expects the input samples to be organized in a two-level hierarchy:

At the first level, input files are grouped by samples (e.g.: patients or biological replicates of an experiment).
At the second level, different datasets belonging to the same sample (e.g., from sample dates) are distinguished.
Inside the 2nd-level directory, the sub-directory `raw_data` holds the sequencing data in FASTQ format (optionally compressed with GZip).
Paired-ended reads need to be in split files with suffixes `_R1` and `_R2`.

```text
📁samples
├───📁patient1
│   └───📁date1
│       └───📁raw_data
│           ├───🧬reads_R1.fastq
│           └───🧬reads_R2.fastq
└───📁patient2
    ├───📁date1
    │   └───📁raw_data
    │       ├───🧬reads_R1.fastq
    │       └───🧬reads_R2.fastq
    └───📁date2
        └───📁raw_data
            ├───🧬reads_R1.fastq
            └───🧬reads_R2.fastq
```

## Preparing a small dataset

In the directory `example_HIV_data` you find a small test dataset that you can run on your workstation or laptop.
The files will have the following structure:

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

## Install V-pipe

After, [having installed V-pipe using the installation tutorial](tutorial_0_install.md), create a new working directory for this analysis:

```bash
cd vp-analysis

# create a new directory and initialise it
mkdir -p work_hiv
cd work_hiv
../V-pipe/init_project.sh

cd ../..
```


## Preparation

Copy the samples directory you created in the step "Preparing a small dataset" to this working directory. (You can display the directory structure with `tree vp-analysis/work_hiv/resources/samples` or `find vp-analysis/work_hiv/resources/samples`.)

```bash
mkdir -p vp-analysis/work_hiv/resources
mv vp-analysis/V-pipe/docs/example_HIV_data/samples vp-analysis/work_hiv/resources/samples
```

Note that:
- by default V-pipe expects its samples in a directory `samples` contained directly in the working directory - i.e. `vp-analysis/work_hiv/sample``
- in this tutorial we put them inside the `resources` subdirectory, and will set the config file accordingly.


### Reference
If you have a reference sequences that you would like to use for read mapping and alignment, then add it to the `resources/reference/ref.fasta` directory. In our case, however, we will use the reference sequence HXB2 already provided by V-Pipe `V-pipe/resources/hiv/HXB2.fasta`.

### Preparing V-pipe's configuration

In the `work_hiv`  directory you can find the file `config.yaml`. This is where the V-Pipe configuation should be specified. See [here](https://github.com/cbg-ethz/V-pipe/tree/master/config#readme) for the documentation of the configuration.
In this tutorial we are building our own configuration therefore `virus_base_config` will remain empty. Since we are working with HIV-1, V-Pipe is providing meta information that will be used for visualisation (metainfo_file and gff_directory).

```bash
cat <<EOT > ./vp-analysis/work_hiv/config.yaml
general:
    virus_base_config: ""
    aligner: bwa
    snv_caller: shorah
    haplotype_reconstruction: haploclique

input:
    reference: "{VPIPE_BASEDIR}/../resources/hiv/HXB2.fasta"
    metainfo_file: "{VPIPE_BASEDIR}/../resources/hiv/metainfo.yaml"
    gff_directory: "{VPIPE_BASEDIR}/../resources/hiv/gffs/"
    # NOTE: this input datadir isn't standard
    datadir: resources/samples/
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
EOT
```

**Note**: A YAML files use spaces as indentation, you can use 2 or 4 spaces for indentation, but **no tab**. There are also [online YAML file validators](https://www.yamllint.com/) that you might want to use if your YAML file is wrongly formatted.

## Running V-pipe


Before running check what will be executed:

```bash
cd vp-analysis/work_hiv/

./vpipe --dryrun

cd ../..
```

As this is your first run of V-pipe, it will also generate the sample collection table. Check `samples.tsv` in your editor.

Note that the samples you have downloaded have reads of length 301 only. V-pipe’s default parameters are optimized for reads of length 250. To adapt to the read length, add a third column in the tab-separated file as follows:

```bash
cat <<EOT > vp-analysis/work_hiv/samples.tsv
CAP217	4390	301
CAP188	4	301
CAP188	30	301
EOT
```

Always check the content of the `samples.tsv` file.

If you did not use the correct directory structure, this file might end up empty or some entries might be missing.
You can safely delete it and re-run with option `--dry-run` to regenerate it.

Finally, we can run the V-pipe analysis (the necessary dependencies will be downloaded and installed in conda environments managed by snakemake):

```bash
cd vp-analysis/work_hiv/

./vpipe -p --cores 2

cd -
```


## Output

The Wiki contains an overview of the output files. The output of the SNV calling step is aggregated in a standard [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) file, located in `results/​{hierarchy}​/variants/SNVs/snvs.vcf`. You can open it with your favorite VCF tools for visualisation or downstream processing. It is also available in a tabular format in `results/​{hierarchy}​/variants/SNVs/snvs.csv`.


## Swapping component

The default configuration uses ShoRAH to call the SNVs and to reconstruct the local (windowed) haplotypes.

Components of the pipeline can be swapped simply by changing the `config.yaml` file. For example to call SNVs using lofreq instead of ShoRAH use

```yaml
general:
  snv_caller: lofreq
```
