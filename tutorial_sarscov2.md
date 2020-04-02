---
layout: default
permalink: /tutorial/sars-cov2/
---

This tutorial shows you the basic of how to interact with V-pipe.

For the purose of this Tutorial, we will work with the `sars-cov2` branch which is [adapted for the SARS-CoV-2 virus]({{ "/sars-cov-2/" | relative_url }}).

> **Organisation of Data** :
>
> V-pipe expects its samples data organized in a two-level hierarchy:
>
> - Input files to be grouped by samples (e.g.: patient samples or biological replicates of an experiment).
> - A second level for distinction of datasets belonging to the same sample (e.g.: sample dates).
> - inside, the directory `raw_data` hold the sequencing output in FASTQ format (optionally compressed with GZip)
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
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# -b for batch (no question asked)
bash Miniconda3-latest-Linux-x86_64.sh -b -p ~/V-test/miniconda3
```

And now we're going to start miniconda and install snakemake from bioconda.

```bash
# mind the dot (=source)
. ~/V-test/miniconda3/bin/activate

conda create -n V-pipe -c bioconda snakemake conda
conda activate V-pipe
# We will let snakemake --use-conda handle
# installation and download of V-pipe dependencies
```
> **Tips:**  To save space, you can install the *-minimal* version of snakemake, without all the GUI dependencies

## Clone V-pipe


```bash
git clone -b sars-cov2 https://github.com/cbg-ethz/V-pipe.git
cd V-pipe
```
> **Versions**: different version come with different defaults
> - *sars-cov2* branch is adapted for SARS-CoV-2
> - the default branch (*master*) is adapted for HIV

## Test V-pipe

Put the `samples` hierarchy that you created before in this `V-pipe` directory. Then verify that it has the [desired structure](#preparing-a-small-dataset) :

```bash
tree samples
```
(or alternatively: `find samples`[^fancyfind])

[^fancyfind]: `find samples | sed -e "s/[^-][^\/]*\// |/g" -e "s/|\([^ ]\)/|-\1/"` is overly long but produces slightly prettier output than `find samples`

Check the parameters of the `vpipe.config` file with you favourite editor (`vim`, `emacs`, `nano`, [butterflies](https://xkcd.com/378/), etc.).
For a SNV and Local (windowed) haplotype reconstruction, you will need at least the following options:

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
V-pipe's defaults parameters are optimized for reads length of 200 ; add the third column in the file:

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


## Larger dataset

Log onto [Euler](https://scicomp.ethz.ch/wiki/Euler) with ETH account.

```bash
ssh login.euler.ethz.ch
```

There is a much larger test data that you can use:

```bash
ls /cluster/work/bewi/members/mpirkl/Sars-Cov-2_data/
ls /cluster/work/bewi/members/mpirkl/Sars-Cov-2_data/PRJNA6*
ls /cluster/work/bewi/members/mpirkl/Sars-Cov-2_data/PRJNA6*/sample*/
```

## Cluster deployment

It is possible to ask snakemake [submit jobs on a cluster](https://snakemake.readthedocs.io/en/stable/executing/cluster-cloud.html#cluster-execution), such as [LSF system used on Euler](https://scicomp.ethz.ch/wiki/Using_the_batch_system)




```bash
# Linux:
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh -b -p $SCRATCH/miniconda3

# Mac OS X:
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
sh Miniconda3-latest-MacOSX-x86_64.sh -b -p $SCRATCH/miniconda3

$SCRATCH/miniconda3/bin/conda install -c bioconda snakemake-minimal
```
> (We can install the *-minimal* version of snakemake: we will probably not run any GUI functionality)

You might get versions conflict at this point, specially if you are re-using your older miniconda version which could be outdated. Consider updating it:

```bash
$SCRATCH/miniconda3/bin/conda update conda
```

And/or if you have a recent enough `conda` executable, you can also consider placing your requirement in a different prefix, like we did in the workstation/laptop section of the tutorial:

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

Check the other [options for running snakemake on clusters](https://snakemake.readthedocs.io/en/stable/executing/cli.html#CLUSTER) if you need more advanced uses.

-----
