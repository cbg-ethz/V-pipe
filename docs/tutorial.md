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

This tutorial shows the basics of how to interact with V-pipe.


## Organizing Data

V-pipe expects the input samples to be organized in a two-level hierarchy:

At the first level, input files grouped by samples (e.g.: patients or biological replicates of an experiment).
A second level for distinction of datasets belonging to the same sample (e.g.: sample dates).
Inside that directory, the sub-directory raw_data holds the sequencing data in FASTQ format (optionally compressed with GZip).
Paired-ended reads need to be in split files with _R1 and _R2 suffixes.

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

You can run the first test on your workstation or a good laptop.

First, you need to prepare the data. As an example, we download the following Multiplexed Illumina MiSeq samples of HIV from SRA: SRR9588829 and SRR9588828. They where taken from a HIV-1 positive patient (with identifier: CAP188) at 4 (SRR9588828) and 30 (SRR9588829) weeks post-infection. Find here the corresponding publicaiton: DOI: 10.1126/scitranslmed.aaw5589


```bash
mkdir -p samples/CAP188/4/
cd samples/CAP188/4/
fastq-dump -O raw_data --split-e  SRR9588828
cd -
```

```bash
mkdir -p samples/CAP188/30/
cd samples/CAP188/30/
fastq-dump -O raw_data --split-e  SRR9588829
cd -
```
Using the `--split-e` option, we download the reads seperatled into forward and reverse reads. Here you can find some more information on `fastq-dump`: https://edwards.flinders.edu.au/fastq-dump/

You then have to rename the files so that they have `_R1` and `_R2` suffixes:

```bash
mv samples/CAP188/4/raw_data/SRR9588828_1.fastq samples/CAP188/4/raw_data/SRR9588828_R1.fastq
mv samples/CAP188/4/raw_data/SRR9588828_2.fastq samples/CAP188/4/raw_data/SRR9588828_R2.fastq

mv samples/CAP188/30/raw_data/SRR9588829_1.fastq samples/CAP188/30/raw_data/SRR9588829_R1.fastq
mv samples/CAP188/30/raw_data/SRR9588829_2.fastq samples/CAP188/30/raw_data/SRR9588829_R2.fastq
```

The downloaded files will have the following structure:

```
samples
└───CAP188
    |───4
    |   └───raw_data
    |       |───SRR9588828_R1.fastq
    |       └───SRR9588828_R2.fastq
    └───30
        └───raw_data
            |───SRR9588829_R1.fastq
            └───SRR9588829_R2.fastq
```


## Install V-pipe

V-pipe uses the [Bioconda](https://bioconda.github.io/) bioinformatics software repository for all its pipeline components. The pipeline itself is written using [Snakemake](https://snakemake.readthedocs.io/en/stable/).

For advanced users: If your are fluent with these tools, you can:

* directly download and install [bioconda](https://bioconda.github.io/user/install.html) and [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda),
* specifiy your V-pipe configuration
* and start using V-pipe with them, using the --use-conda to [automatically download and install](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management) any further pipeline dependencies.
* please refer to the documentation for additional instructions.

The present tutorial will show simplified commands that automate much of this process.

To deploy V-pipe, you can use the installation script with the following parameters:

```bash
curl -O 'https://raw.githubusercontent.com/cbg-ethz/V-pipe/master/utils/quick_install.sh'
bash quick_install.sh -p testing -w work
```
If you get `zsh: permission denied: ./quick_install.sh`, then run `chmod +x quick_install.sh` this gives the necessary permissions.

* using `-p` specifies the subdirectory where to download and install snakemake and V-pipe
* using `-w` will create a working directory and populate it. It will copy over the references and the default `config/config.yaml`, and create a handy `vpipe` short-cut script to invoke `snakemake`.

Tip: To create and populate other new working directories, you can call init_project.sh from within the new directory:

```console
mkdir -p working_2
cd working_2
../V-pipe/init_project.sh
```

## Preparation

Copy the samples directory you created in the step Preparing a small dataset to this working directory. You can display the directory structure with tree samples or find samples.

```bash
mv ./samples ./testing/work/resources/
```

### Reference
If reference sequences is available add it to the `resources/reference/ref.fasta` directory.
In our case, we download the referecne sequecne: https://www.ncbi.nlm.nih.gov/nuccore/MN097551.1

### Preparing V-pipe's configuration

In the `work`  directory you can find the file `config.yaml`. This is where the V-Pipe configuation should be specified. See [here] (https://github.com/cbg-ethz/V-pipe/tree/master/config#readme) for the documentation of the configuation. In this tutorial we are building our own configuarion therefore `virus_base_config` will remain empty.

```bash
general:
    virus_base_config: ''
    aligner: "bwa"
    snv_caller: "lofreq"
    haplotype_reconstruction: "haploclique"

input:
    reference: "{VPIPE_BASEDIR}/../../work/resources/reference/ref.fasta"
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
    visualization: false
    QA: false
    diversity: false
```

Tip: A YAML file use spaces as indentation, you can use 2 or 4 spaces for indentation, but no tab. There is also online yaml file validator that you might want to check if your yaml file is wrongly formatted.

## Running V-pipe


Before running heck what will be executed:

```bash
cd ./testing/work/
./vpipe --dryrun
```

As it is your first run of V-pipe, this will also generate the sample collection table. Check `samples.tsv` in your editor.

Note that the demo files you downloaded have reads of length 150 only. V-pipe’s default parameters are optimized for reads of length 250; add the third column in the tab-separated file:

```bash
cat <<EOT > ./testing/work/samples.tsv
SRR10903401	20200102	301
SRR10903402	20200102	301
EOT
```

Tip: Always check the content of the `samples.tsv` file.

If you didn’t use the correct structure, this file might end up empty or some entries might be missing.
You can safely delete it and re-run the `--dryrun` to regenerate it.

Run the V-pipe analysis (the necessary dependencies will be downloaded and installed in conda environments managed by snakemake):

```bash
cd ./testing/work/
./vpipe -p --cores 2
```


## Output

The Wiki contains an overview of the output files. The output of the SNV calling is aggregated in a standard [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) file, located in `samples/​{hierarchy}​/variants/SNVs/snvs.vcf`, you can open it with your favorite VCF tools for visualisation or downstream processing. It is also available in a tabular format in `samples/​{hierarchy}​/variants/SNVs/snvs.csv`.

### Expected output

The small dataset that we used in this tutorial section has been analyzed by [doi:10.1093/nsr/nwaa036](https://doi.org/10.1093/nsr/nwaa036). The results of the original analysis (using bwa, samtools mpileup, and bcftools) are displayed in Table 2 in the article:

Using either the VCF or CSV files, compare with the results given out by V-pipe (with bwa and ShoRAH).

* For positions 19164 and 24323 of SRR10903401 and position 11563 of SRR10903402, we expect to see similar results in V-pipe.
* For the remaining positions (1821, 26314 and 26590 of SRR10903401), we expect that ShoRAH will consider the variants of poor quality and reject them because there is very little support ( <= than 5 reads supporting the alt).


## Swapping component

The default configuration uses ShoRAH to call the SNVs and to reconstruct the local (windowed) haplotypes.

Components can be swapped simply by changing the `config.yaml` file. For example to call SNVs using lofreq:

```yaml
general:
  snv_caller: lofreq
```
