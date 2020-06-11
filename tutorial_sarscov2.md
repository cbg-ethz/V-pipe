---
layout: default
permalink: /tutorial/sars-cov2/
---

This tutorial shows the basics of how to interact with V-pipe.

For the purpose of this Tutorial, we will work with the `sars-cov2` branch which is [adapted for the SARS-CoV-2 virus]({{ "/sars-cov-2/" | relative_url }}).

> **Organizing Data**:
>
> V-pipe expects [the input samples to be organized in a two-level hierarchy](https://github.com/cbg-ethz/V-pipe/wiki/getting-started#input-files):
>
> - At the first level, input files grouped by samples (e.g.: patients or biological replicates of an experiment).
> - A second level for distinction of datasets belonging to the same sample (e.g.: sample dates).
> - Inside that directory, the sub-directory `raw_data` holds the sequencing data in FASTQ format (optionally compressed with GZip).
> - Paired-ended reads need to be in split files with `_R1` and `_R2` suffixes.

*[FASTQ]: [represents DNA sequencer reads along with quality scores](https://en.wikipedia.org/wiki/FASTQ_format)

## Preparing a small dataset

You can run the first test on your workstation or a good laptop.

First, you need to prepare the data:

 - For that test, you need to download the following runs from SRA:
   [SRR10903401](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR10903401) and
   [SRR10903402](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR10903402)
   
   (hint: click on the _data_ tab to find the download links).
 - Then, organize the downloaded files to have the following structure
   (you need to rename the files so that they have `_R1` and `_R2` suffixes):


*[SRA]: [sequence read archive](https://trace.ncbi.nlm.nih.gov/Traces/sra/)

```
samples
├── SRR10903401
│   └── 20200102
│       └── raw_data
│           ├── wuhan2_R1.fastq
│           └── wuhan2_R2.fastq
└── SRR10903402
    └── 20200102
        └── raw_data
            ├── wuhan1_R1.fastq
            └── wuhan1_R2.fastq
```

## Install V-pipe

V-pipe uses the [Bioconda](https://bioconda.github.io/)[^bioconda] bioinformatics software repository for all its pipeline components.
The pipeline itself is written using [snakemake](https://snakemake.readthedocs.io/)[^Snakemake].

[^bioconda]: Grüning, Björn, Ryan Dale, Andreas Sjödin, Brad A. Chapman, Jillian Rowe, Christopher H. Tomkins-Tinch, Renan Valieris, the Bioconda Team, and Johannes Köster. 2018. “Bioconda: Sustainable and Comprehensive Software Distribution for the Life Sciences”. Nature Methods, 2018 doi:[10.1038/s41592-018-0046-7](https://doi.org/10.1038/s41592-018-0046-7).

[^Snakemake]: Johannes Köster and Sven Rahmann. [Snakemake – a scalable bioinformatics workflow engine](https://academic.oup.com/bioinformatics/article/28/19/2520/290322). *Bioinformatics*, 28(19):2520–2522, 2012. doi:[10.1093/bioinformatics/bts480](https://doi.org/10.1093/bioinformatics/bts480)

> **For advanced users:** If your are fluent with these tools, you can:
> - directly download and install
>   [bioconda](https://bioconda.github.io/user/install.html) and
>   [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda),
> - clone the [sars-cov2 branch of V-pipe](https://github.com/cbg-ethz/V-pipe/tree/sars-cov2)
> - and start using V-pipe with them, using the `--use-conda` to
>   [automatically download and install](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management)
>   any [further pipeline dependencies]({{ "/pipeline/" | relative_url }}).
> - please refer to [the documentation](https://github.com/cbg-ethz/V-pipe/wiki/getting-started#running-v-pipe) for additional instructions.
>
> The present tutorial will show simplified commands that automate much of this process.


To deploy V-pipe, you can use the installation script with the following parameters:

```bash
curl -O 'https://raw.githubusercontent.com/cbg-ethz/V-pipe/master/utils/quick_install.sh'
bash quick_install.sh -b sars-cov2 -p testing -w work
cd ./testing/work/
```

 - using `-p` specifies the subdirectory where to download and install snakemake and V-pipe
 - using `-b` specifies which branch of V-pipe to use

   > **Versions**: different versions come with different defaults.
   > - `sars-cov2` branch is adapted for SARS-CoV-2
   > - the default branch (`master`) is adapted for HIV
   >
   > As the project matures, [release tar-balls frozen at specific versions](https://github.com/cbg-ethz/V-pipe/releases) will progressively be made available, see the [Cluster example below](#cluster-deployment).
   >
 - using `-w` will create a working directory and populate it. It will copy over the references and the default `vpipe.config`, and create a handy `vpipe` short-cut script to invoke `snakemake`.

   > **Tips**:
   > To create and populate other new working directories, you can call `init_project.sh` from within the new directory :
   ```bash
   mkdir -p working_2
   cd working_2
   ../V-pipe/init_project.sh
   ```

## Running V-pipe


Copy the `samples` directory you created in the step [Preparing a small dataset](#preparing-a-small-dataset) to this `working` directory. 
You can display the directory structure with `tree samples` or `find samples`[^fancyfind].

[^fancyfind]: `find samples | sed -e "s/[^-][^\/]*\// |/g" -e "s/|\([^ ]\)/|-\1/"` is overly long but produces slightly prettier output than `find samples`

Check the parameters of the `vpipe.config` file with you favorite editor (`vim`, `emacs`, `nano`, [butterflies](https://xkcd.com/378/), etc.).
For SNVs and Local (windowed) haplotype reconstruction, you will need at least the following options:

*[SNVs]: Single Nucleotide Variants

```
[input]
reference = references/NC_045512.2.fasta

[output]
snv = True
local = True
global = False

[general]
aligner = bwa
```

Check what will be executed:
```bash
./vpipe --dryrun
```

As it is your first run of V-pipe, this will also generate the sample collection table.
Check `samples.tsv` in your editor.

Note that the demo files you downloaded have reads of length 150 only.
V-pipe's default parameters are optimized for reads of length 250 ; add the third column in the tab-separated file:

```
SRR10903401	20200102	150
SRR10903402	20200102	150
```

> **Tips:** Always check the content of the `samples.tsv` file.
>  - If you didn't use the [correct structure](#preparing-a-small-dataset), this file might end up empty or some entries might be missing.
>  - You can safely delete it and re-run the `--dryrun` to regenerate it.

Run the V-pipe analysis (the necessary dependencies will be downloaded and installed in conda environments managed by snakemake):

```bash
./vpipe -p --cores 2
```

## Output

The Wiki contains an overview of the [output files](https://github.com/cbg-ethz/V-pipe/wiki/output).
The output of the SNV calling is aggregated in a standard [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) file, located in
`samples/`*\{hierarchy\}*`/variants/SNVs/snvs.vcf`, you can open it with your favorite VCF tools for visualisation or downstream processing.
It is also available in a tabular format in `samples/`*\{hierarchy\}*`/variants/SNVs/snvs.csv`.

*[VCF]: Variant Call Format


> **Note:** The visualization and reporting features are still being continuously updated.

### Expected output

The small dataset that we used in this tutorial section has been analyzed by [doi:10.1093/nsr/nwaa036](https://doi.org/10.1093/nsr/nwaa036).
The results of the original analysis (using bwa, samtools mpileup, and bcftools) are displayed in Table 2 in the article:

<div class="table-wrapper" markdown="block">

|Accession number|Genomic position|Ref allele|Alt allele|Ref reads|Alt reads|Location_date  |GISAID ID     |
|:---------------+---------------:+:--------:+:--------:+--------:+--------:+:--------------+:-------------|
|SRR10903401     |            1821|     G    |     A    |       52|        5|WH_2020/01/02.a|EPI_ISL_406716|
|SRR10903401     |           19164|     C    |     T    |       40|       12|WH_2020/01/02.a|EPI_ISL_406716|
|SRR10903401     |           24323|     A    |     C    |      102|       67|WH_2020/01/02.a|EPI_ISL_406716|
|SRR10903401     |           26314|     G    |     A    |       15|        2|WH_2020/01/02.a|EPI_ISL_406716|
|SRR10903401     |           26590|     T    |     C    |       10|        2|WH_2020/01/02.a|EPI_ISL_406716|
|SRR10903402     |           11563|     C    |     T    |      164|       26|WH_2020/01/02.b|EPI_ISL_406717|

</div>

Using either the VCF or CSV files, compare with the results given out by V-pipe (with bwa and ShoRAH).

- For positions 19164 and 24323 of SRR10903401 and position 11563 of SRR10903402,
  we expect to see similar results in V-pipe.
- For the remaining positions (1821, 26314 and 26590 of SRR10903401),
  we expect that ShoRAH will consider the variants of poor quality and reject them
  because there is very little support ( \<=  than 5 reads supporting the alt).

## Swapping component

The default configuration uses [ShoRAH](https://cbg-ethz.github.io/shorah/)
to call the SNVs and to reconstruct the local (windowed) haplotypes.

Components can be swapped simply by changing the `vpipe.config` file.
For example to call SNVs using `lofreq`:

```
[output]
snv = True
local = False

[general]
snv_caller=lofreq
```


## Cluster deployment

It is possible to ask snakemake to [submit jobs on a cluster](https://snakemake.readthedocs.io/en/stable/executing/cluster-cloud.html#cluster-execution) using the batch submission command-line interface of your cluster.

[Platform LSF by IBM](https://www.ibm.com/support/knowledgecenter/SSETD4_9.1.2/lsf_command_ref/bsub.1.html) is one of the popular systems you might find (Others include [SLURM](https://slurm.schedmd.com/sbatch.html), [Grid Engine](https://en.wikipedia.org/wiki/Oracle_Grid_Engine)).

*[LSF]: Load Sharing Facility

To deploy on the cluster:

```bash
wget 'https://raw.githubusercontent.com/cbg-ethz/V-pipe/master/utils/quick_install.sh'
bash quick_install.sh -b sars-cov2 -p $SCRATCH -w working
cd $SCRATCH/working/
```

 - using `-p` help us store V-pipe on some large-storage share, e.g.: scratch suffice for this tutorial
 
> **Tips:** As V-pipe for SARS-CoV-2 matures, it will be possible to download [snapshots frozen at specific version](https://github.com/cbg-ethz/V-pipe/releases).
> This enables more reproducible results. To specify a release use the `-r` option :
```bash
wget 'https://raw.githubusercontent.com/cbg-ethz/V-pipe/master/utils/quick_install.sh'
bash quick_install.sh -r sars-cov2-snapshot-20200406 -p $SCRATCH -w working
cd $SCRATCH/working/
```
> this will download the tarball [sars-cov2-snapshot-20200406.tar.gz](https://github.com/cbg-ethz/V-pipe/archive/sars-cov2-snapshot-20200406.tar.gz) and uncompress it into a directory called `V-pipe-sars-cov2-snapshot-20200406`

### Running V-pipe on the cluster

In the `working` directory, create a `samples` sub-directory and populate it. Check its structure with `tree` or `find`. Perform the necessary adjustments to `vpipe.config`.

To run V-pipe on a cluster :

> **Tips:** There are [snakemake parameters for conda](https://snakemake.readthedocs.io/en/stable/executing/cli.html#CONDA)
that can help management of dependencies:
> - using `--conda-create-envs-only` enables to download the dependencies only without running the pipeline itself.
> - using `--conda-prefix {DIR}`  stores the conda environments of dependencies in a common directory (thus possible to share re-use between multiple instances of V-pipe).
>
> When using V-pipe in production environments, plan the `-p` prefix, `-w` working and `--conda-prefix` environments directories according to the cluster quotas and time limits

```bash
# Download everything in advance
./vpipe --conda-prefix $SCRATCH/snake-envs --cores 1 --conda-create-envs-only

# Cluster LSF submitting
./vpipe --conda-prefix $SCRATCH/snake-envs -p --cluster 'bsub' --jobs 2

# Using bsub on the master job too, instead of running it on the login node
bsub ./vpipe --conda-prefix $SCRATCH/snake-envs -p --cluster 'bsub' --jobs 2

# Alternative for running everything from a single interactive SSH node
bsub -I <<<"./vpipe --conda-prefix $SCRATCH/snake-envs -p --cores 2"
```
> **Tips:** See the V-pipe documentation for [more cluster commands](https://github.com/cbg-ethz/V-pipe/wiki/advanced#running-v-pipe-on-a-lsf-cluster).


Check the other [options for running snakemake on clusters](https://snakemake.readthedocs.io/en/stable/executing/cli.html#CLUSTER) if you need more advanced uses.

-----