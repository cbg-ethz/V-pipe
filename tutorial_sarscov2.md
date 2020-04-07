---
layout: default
permalink: /tutorial/sars-cov2/
---

This tutorial shows you the basic of how to interact with V-pipe.

For the purpose of this Tutorial, we will work with the `sars-cov2` branch which is [adapted for the SARS-CoV-2 virus]({{ "/sars-cov-2/" | relative_url }}).

> **Organizing Data** :
>
> V-pipe expects its [samples data organized in a two-level hierarchy](https://github.com/cbg-ethz/V-pipe/wiki/getting-started#input-files):
>
> - Input files to be grouped by samples (e.g.: patient samples or biological replicates of an experiment).
> - A second level for distinction of datasets belonging to the same sample (e.g.: sample dates).
> - Inside, the directory `raw_data` hold the sequencing output in FASTQ format (optionally compressed with GZip)
> - When in split files, paired-ends reads need to have `_R1` and `_R2` suffixes in their name.

*[FASTQ]: [represents DNA sequencer reads along with quality scores](https://en.wikipedia.org/wiki/FASTQ_format)

## Preparing a small dataset

You can run the first test on your workstation or a good laptop.

First you'll need to prepare the data:

 - For that test, you need to download from SRA the following runs:
   [SRR10903401](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR10903401) and
   [SRR10903402](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR10903402)
   
   (hint: click on the _data_ tab)
 - Then arrange the downloaded file to have the following structure
   (you need to rename the files so that the have `_R1` and `_R2` suffixes) :


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

## Miniconda

V-pipe uses the [Bioconda](https://bioconda.github.io/)[^bioconda] bioinformatics software repository for all component of its pipeline.
The pipeline itself is written using [snakemake](https://snakemake.readthedocs.io/)[^Snakemake].

[^bioconda]: Grüning, Björn, Ryan Dale, Andreas Sjödin, Brad A. Chapman, Jillian Rowe, Christopher H. Tomkins-Tinch, Renan Valieris, the Bioconda Team, and Johannes Köster. 2018. “Bioconda: Sustainable and Comprehensive Software Distribution for the Life Sciences”. Nature Methods, 2018 doi:[10.1038/s41592-018-0046-7](https://doi.org/10.1038/s41592-018-0046-7).

[^Snakemake]: Johannes Köster and Sven Rahmann. [Snakemake – a scalable bioinformatics workflow engine](https://academic.oup.com/bioinformatics/article/28/19/2520/290322). *Bioinformatics*, 28(19):2520–2522, 2012. doi:[10.1093/bioinformatics/bts480](https://doi.org/10.1093/bioinformatics/bts480)

We're going to set up out environment:

```bash
# directory for this tutorial
mkdir V-test
cd V-test

# Download Miniconda3

# Linux:
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh -b -p ~/V-test/miniconda3

# Mac OS X:
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
sh Miniconda3-latest-MacOSX-x86_64.sh -b -p ~/V-test/miniconda3

# -b for batch (no question asked)
```

And now we're going to start miniconda and install snakemake from bioconda.

```bash
# mind the dot (=source)
. ~/V-test/miniconda3/bin/activate

conda create -n V-pipe -c conda-forge -c bioconda snakemake conda
conda activate V-pipe
# We will let snakemake --use-conda handle
# installation and download of V-pipe dependencies
```
> **Tips:**  To save space, you can install the *-minimal* version of *snakemake*, without all the GUI dependencies.
> If you don't have any development tools on this machine, you can also install *git* and your favorite text editor

## Clone V-pipe


```bash
git clone https://github.com/cbg-ethz/V-pipe.git
cd V-pipe
git checkout sars-cov2
```
> **Versions**: different version come with different defaults.
> - *sars-cov2* branch is adapted for SARS-CoV-2
> - the default branch (*master*) is adapted for HIV
>
> As the project matures, [release tar-balls frozen at specific versions](https://github.com/cbg-ethz/V-pipe/releases) will progressively be made available, see the [Cluster example below](#Cluster-deployment).

## Test V-pipe

Put the `samples` hierarchy that you created before in this `V-pipe` directory. Then verify that it has the [desired structure](#preparing-a-small-dataset) :

```bash
tree samples
```
(or alternatively: `find samples`[^fancyfind])

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
snakemake -s vpipe.snake --dryrun -p --cores 2
```

As it is your first run of V-pipe, this will also generate the sample collection table.
Check `samples.tsv` in your editor.

Note that the demo files you downloaded have only have reads of length 150.
V-pipe's defaults parameters are optimized for reads length of 250 ; add the third column in the file:

```
SRR10903401	20200102	150
SRR10903402	20200102	150
```

> **Tips:** Always check the content of the `samples.tsv` file.
>  - If you didn't use the [correct structure](#preparing-a-small-dataset) this file might end up empty or some entries might be missing.
>  - You can safely delete it and re-run the `--dryrun` to regenerate it

Run the V-pipe analysis (the necessary dependencies will be downloaded and installed in conda environments managed by snakemake):

```bash
snakemake -s vpipe.snake --use-conda -p --cores 2
```

## Output

<< This section still needs polishing before final release on https://cbg-ethz.github.io/V-pipe/ >>

The Wiki contains an overview of the [output files](https://github.com/cbg-ethz/V-pipe/wiki/output).
The output of the SNVs calling is aggregated in a standard [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) file, located in
`samples/`*\{hierarchy\}*`/variants/SNVs/snvs.vcf`, you can open it with your favorite VCF tools for aggregated or downstream processing.

*[VCF]: Variant Call Format

There is also a visual report generated in `samples/`*\{hierarchy\}*`/visualization/index.html` that you can simply open with your favorite web browser (Firefox, Chrome). You will get an output similar to that:

![SARS-CoV-2 report visualisation]({{ "/img/visualisation-sars-cov2.png" | relative_url }})

You can use the mouse selection to zoom in region and double clicks to zoom out.

### Expected output

The small dataset that we used in this first part has been analyzed by [doi:10.1093/nsr/nwaa036](https://doi.org/10.1093/nsr/nwaa036).
We find their results (analyzed with bwa, samtools mpileup and bcftools) in the table 2 of the article:

|Accession number|Genomic position|Ref allele|Alt allele|Ref reads|Alt reads|Location_date  |GISAID ID     |
|:---------------+---------------:+:--------:+:--------:+--------:+--------:+:--------------+:-------------|
|SRR10903401     |            1821|     G    |     A    |       52|        5|WH_2020/01/02.a|EPI_ISL_406716|
|SRR10903401     |           19164|     C    |     T    |       40|       12|WH_2020/01/02.a|EPI_ISL_406716|
|SRR10903401     |           24323|     A    |     C    |      102|       67|WH_2020/01/02.a|EPI_ISL_406716|
|SRR10903401     |           26314|     G    |     A    |       15|        2|WH_2020/01/02.a|EPI_ISL_406716|
|SRR10903401     |           26590|     T    |     C    |       10|        2|WH_2020/01/02.a|EPI_ISL_406716|
|SRR10903402     |           11563|     C    |     T    |      164|       26|WH_2020/01/02.b|EPI_ISL_406717|

Using the reports and zoom function, try to compare with the results given out by V-pipe (with bwa and ShoRAH).

- For positions 19164 and 24323 of SRR10903401 and position 11563 of SRR10903402,
  we expect to see similar results in V-pipe.
- For the remaining position (1821 of SRR10903401 and 26314, 26590 SRR10903402),
  as there is very little support ( \<=  than 5 reads supporting the alt)
  we expect that ShoRAH will consider the variants of poor quality and reject them.

## Swapping component

The default configuration uses [ShoRAH](https://cbg-ethz.github.io/shorah/)
to call the SNVs and compute the Local (windowed) haplotypes.

You can swap component simply by changing the `vpipe.config` file.
For example to compute SNVs using `lofreq`:

```
[output]
snv = True
local = False

[general]
snv_caller=lofreq
```


## Cluster deployment

It is possible to ask snakemake [submit jobs on a cluster](https://snakemake.readthedocs.io/en/stable/executing/cluster-cloud.html#cluster-execution) using the batch submission command-line interface of your cluster.

[Platform LSF by IBM](https://www.ibm.com/support/knowledgecenter/SSETD4_9.1.2/lsf_command_ref/bsub.1.html) is one of the popular system you might find. (Others include [SLURM](https://slurm.schedmd.com/sbatch.html), [Grid Engine](https://en.wikipedia.org/wiki/Oracle_Grid_Engine) )

*[LSF]: Load Sharing Facility


```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh -b -p $SCRATCH/miniconda3

$SCRATCH/miniconda3/bin/conda install -c bioconda snakemake-minimal
```
> (We can install the *-minimal* version of snakemake: we will probably not run any GUI functionality)

Specially if you are re-using your older miniconda version which could be outdated, you might get versions conflict at this point. Consider updating it:

```bash
$SCRATCH/miniconda3/bin/conda update conda
```

Recent miniconda environment allow you to install conda updates in different prefixes/conda environments.
You can create a separate conda environment prefix with everything you need to start V-pipe (in the same way we set it up in the earlier example) :

```bash
$SCRATCH/miniconda3/bin/conda create -p $SCRATCH/V-pipe_conda -c bioconda conda snakemake
$SCRATCH/V-pipe_conda/bin/snakemake --version
```

Let's fetch V-pipe:

```bash
cd $SCRATCH
git clone -b sars-cov2 https://github.com/cbg-ethz/V-pipe.git
cd V-pipe
```
> **Tips:** as V-pipe for SARS-CoV-2 matures, it will be possible to download snapshots frozen at specific version, for more reproducible results.

```bash
wget https://github.com/cbg-ethz/V-pipe/archive/sars-cov2-snapshot-20200406.tar.gz
tar xvzf sars-cov2-snapshot-20200406.tar.gz
cd V-pipe-sars-cov2-snapshot-20200406/
```

> **Tips:** There are [snakemake parameters for conda](https://snakemake.readthedocs.io/en/stable/executing/cli.html#CONDA)
that can help manage:
> - using `--create-envs-only`, it's possible to only download the dependencies, but not run the pipeline itself
> - using `--conda-prefix {DIR}`  will store the conda environments of dependencies in that directory (thus possible to share re-use between multiple instances of V-pipe)

```bash
. $SCRATCH/miniconda3/bin/activate

# download everything in advance
$SCRATCH/miniconda3/bin/snakemake -s vpipe.snake --use-conda --conda-prefix $SCRATCH/snake-envs --cores all --create-envs-only

# cluster LSF dispatching
$SCRATCH/miniconda3/bin/snakemake -s vpipe.snake --use-conda --conda-prefix $SCRATCH/snake-envs -p --cluster 'bsub' --jobs 2

# alternative for running everything from a single interactive SSH node
bsub -I <<<"$SCRATCH/miniconda3/bin/snakemake -s vpipe.snake --use-conda --conda-prefix $SCRATCH/snake-envs -p --cores 2 --jobs 2"
```
> **Tips:** see the documentation for [more cluster commands](https://github.com/cbg-ethz/V-pipe/wiki/advanced#running-v-pipe-on-a-lsf-cluster).


Check the other [options for running snakemake on clusters](https://snakemake.readthedocs.io/en/stable/executing/cli.html#CLUSTER) if you need more advanced uses.

-----
