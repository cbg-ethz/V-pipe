---
layout: page
title: Tutorial
permalink: /tutorial/
---

# TODO

* NEXT DATES : BC2 2019
* Figures
* Modernize image to OVA
 * The download links are provided via CISD File EXchanger (CIFEX).
* alternative: container

-----


In this tutorial, we will learn how to use V-pipe for analyzing NGS data of short viral genomes. The workflow
is designed for both genomics research, as well as clinical diagnostics of RNA and DNA viruses. Among other
features, V-pipe assesses viral diversity on three genomic scales:

  - (i) SNV: Frequencies and positions of singlenucleotide variants that differ from the population
  - (ii) Local: co-occurrence of SNVs in regions that areas long as the average read
  - (iii) Global: haplotypes of larger segments of viral genomes.

Moreover, the pipeline can be applied for detecting flowcell cross contamination from other sequencing runs.
V-pipe is available in the [CBG-ETH github repository](https://github.com/cbg-ethz/V-pipe), and we encourage
users to use [this platform for discussing](/contact/) potential issues.

V-Pipe uses Snakemake[^Snakemake], a robust workflow management system written in Python. Therefore, we start
by introducing important concepts and features of Snakemake. Then, we execute the pipeline on downsampled
data[^1] and learn how to adjust the default settings.

However, first of all, we need to set up the environment for running Vpipe. For the tutorial, we provide
a virtual machine (VM) for use with [Oracle VirtualBox](https://www.virtualbox.org/).  This approach is suitable
for Linux, Mac and Windows. While this tutorial assumes that you use the VM, instructions for installing Vpipe
directly on a UNIX-like system are provided in [the appendix](#a-alternative-installation-for-unix-like-operating-systems).


[^1]: This is done in order to reduce waiting times

[^Snakemake]: Johannes Köster and Sven Rahmann. [Snakemake – a scalable bioinformatics workflow engine](https://academic.oup.com/bioinformatics/article/28/19/2520/290322). *Bioinformatics*, 28(19):2520–2522, 2012. doi:[10.1093/bioinformatics/bts480](https://doi.org/10.1093/bioinformatics/bts480)



## 1. Recommended set up for this tutorial: Virtual Machine

1. Download and install [Oracle VirtualBox](https://www.virtualbox.org/)
   (see [documentation](https://www.virtualbox.org/wiki/Downloads)
   ).
2. Download the OVA file `V-pipe.ova`
3. Download and uncompress the test data and scripts for the tutorial.
   Save them to a known location, we recommend using `/home/<user>/vboxShared`.
4. Open the VirtualBox manager. Go to the `File` menu and select `Import Appliance... Ctrl+I`.
   Click on the folder icon (*Choose a virtual appliance file to import...*) next to the input field.
   Then, select the location where you have saved the OVA file.
   Select the `V-pipe.ova` and click on Open.
   You should now see `V-pipe` as a registered VM in the left panel of the VirtualBox Manager window.
5. In order to prevent the image from becoming unmanageably large, we configure a shared folder.
   **Before starting the VM**, you need to change the `host path` to the directory where you have saved
   the tutorial scripts and data. To do so, follow the following steps:

   - (i) In the VirtualBox Manager, select the newly-added VM.
   - (ii) Select `Settings` and then `Shared Folders`.
   - (iii) Under the folders list, double click on the pre-configured shared folder and change the folder
     path to the path where you have saved and uncompressed the tutorial folder.

   Within the VM, the contents of this folder will be accessible via `/media/sf_vboxShared`

6. To test the VM, you only need to click on Start. Login as user `user` using password `vpipe2017`.

**NOTE:** Remember to shut down the virtual machine by clicking on switch-off icon located on the bottom-
right corner. Simply quitting VirtualBox without shutting down the VM may cause a lenghtly filesystem check 
on the next start.

*[OVA]: Open Virtualization Appliance
*[VM]: Virtual Machine



## 2. Snakemake in a nutshell

Since Snakemake is inspired by the GNU make utility, the Vpipe workflow consists of so-called rules. Rules
define how to obtain output files from input files. Rules can depend on other rules, which enables dependencies
to be managed transparently by Snakemake (see [documentation](https://snakemake.readthedocs.io/en/stable/) 
for further information)


### 2.1. "Hello V-pipe" – create your first workflow

**Exercise** Let’s start with a simple exercise. Open your favorite text editor (e.g.: `vim`, `vi`, or `nano`)
and type the following:

```  python
rule all:
	input: "hello_vpipe.txt"

rule hello:
	output: "hello_vpipe.txt"
	shell:
		"""
		echo Hello V-pipe! > {output}
		"""
```

Save the workflow as `hello_vpipe.snake`. In order to run your first workflow, open a terminal and change
to the directory to which the aforementioned file was saved. Then, type:

``` shell
snakemake --snakefile hello_vpipe.snake 
```

By default, Snakemake searches for a file named `Snakefile`. The option `--snakefile`, or its equivalent short
form `-s`, can be used to specify a different *Snakefile* file name, as shown in this example.

If your workflow ran successfully, you should see an output message similar to the following:

```
Provided cores: 1
Rules claiming more threads will be scaled down.
Job counts:
	count 	jobs
	1	all
	1	hello
	2
	
rule hello:
	output: hello_vpipe.txt
	jobid: 1

Finished job 1.
1 of 2 steps (50%) done

localrule all:
	input: hello_vpipe.txt
	jobid: 0

Finished job 0.
2 of 2 steps (100%) done
```

From this output, we can verify that two jobs, corresponding to rules `all` and `hello` , were successfully
submitted and executed. In addition to the output message, a file named `hello_vpipe.txt` should have
been created containing the following text:


    Hello V-pipe!


This example consists of two rules, namely rules `all` and `hello`. In general, a rule is defined by a name
and some directives. Here, we use the directives `input`, `output`, and `shell`. Other directives include `run`,
`script`, `params`, `benchmark`, `threads` and `conda`. Once we come to use these directives, more details will be
provided.

For every rule, directives `shell`, `run`, or `script` indicate how to obtain the output(s) from inputs, if the latter
are specified. As can be deduced from the names, `shell`, `run`, and `script` directives are followed by shell commands, 
Python codes, or scripts (e.g.: R or Python scripts), respectively. For example, `hello_vpipe.snake` 
can be rewritten using Python code as follows:

``` python
rule all:
	input:"hello_vpipe.txt"

rule hello:
	output: "hello_vpipe.txt"
	run:
		with open(output[0], 'w') as fout:
		fout.write("Hello V pipe!")
```

Note that variables are accessed differently when using `shell` and `run` directives. In the former case, curly
braces are employed, while in the latter, variables are accessed as Python lists.

Snakemake determines the order in which rules are executed in a top-down manner.  In this example,
Snakemake will start by searching how to generate the target file `hello_vpipe.txt`, finding out that rule
`hello` needs to be executed first, followed by rule `all` (see Figure 1). In particular, the rule named
`all`[^2] is an special rule. Being the first rule of the Snakefile, it defines what the target files are
(i.e.: the final output of your pipeline).

Try to run `hello_vpipe.snake` once more. You should see the following message:


    Nothing to be done.


This is because, the target file has been already produced.

For more information about Snakemake, please refer to the [documentation](https://snakemake.readthedocs.io/en/stable/) 
and [tutorial](http://slides.com/johanneskoester/snakemake-tutorial-2016#/).


[^2]: The rule `all` can take any name, but is a common practice to name it as such.


### 2.2. Wildcards for rule generalization

It is very desirable, if not a must, to be able to apply the pipeline to different datasets. This is possible using wildcards. 
A wildcard, similar to the asterisk in Bash, represents a slot in which any value could fit. Using wildcards, hard-coded variable
names for input and output files can be avoided, and rules are g generalized to new datasets.

**Exercise** We will use wildcards to create a pipeline for reconstructing viral haplotypes from NGS data
using software Haploclique[^Haploclique]. In a text editor, type the following:

``` python
rule haploclique:
	input:
		"{sample}_aln.bam"
	output:
		"{sample}_quasispecies.fasta"
	params:
		CLIQUE_SIZE_LIMIT = 3,
		MAX_NUM_CLIQUES = 10000,
		HAPLOCLIQUE = "/path/to/haploclique"
	shell:
		"""
		{params.HAPLOCLIQUE} --no_singletons --no_prob0 --limit_clique_size={params.CLIQUE_SIZE_LIMIT} --max_cliques={params.MAX_NUM_CLIQUES} --bam {input} {wildcards.sample}_quasispecies
		"""
```

[^Haploclique]: Armin Töpfer, Tobias Marschall, Rowena A. Bull, Fabio Luciani, Alexander Schönhuth, and Niko Beerenwinkel. [Viral quasispecies assembly via maximal clique enumeration.](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003515) *PLOS Computational Biology*, 10(3):1–10,032014. doi:[10.1371/journal.pcbi.1003515](https://doi.org/10.1371/journal.pcbi.1003515)

For this exercise, the path to the software Haploclique needs to be specified. However, in the VM, Haploclique
is already in the `PATH`, such that this parameter can be omitted.

Here, we use the wildcard `sample` to represent many different inputs. Note that wildcards are written in curly
braces in the `input`, `output`, and `params` directives, while they are accessed using the keyword `wildcards`
in the `shell` directive. In this exercise, we use a new directive, namely `params`. Specifying parameters using
the `params` directive is particularly useful when running the pipeline in a cluster environment (See [section 7](7-running-v-pipe-on-a-lsf-cluster)
), as well for encapsulating user-configurable options.



The test datasets for this exercise are located in the `eg_haploclique` folder. Open a terminal, change into
the directory to which you have saved the scripts and the test datasets (e.g. inside the VM, 
`/media/sf_vboxShared/eg_haploclique`). Then, save the pipeline as `haploclique.snake`. In order to execute it,
type on the terminal:

``` shell
snakemake -s haploclique.snake 2VM-sim_quasispecies.fasta -p
```

In this case, it is required to indicate the target file, e.g.: `2VM-sim_quasispecies.fasta`, in the execution
command. In addition, we use option `-p` for printing the shell commands that are executed. Below, we show
a snippet of the output generated while running the script:

```
rule haploclique:
    input: 2VM-sim_aln.bam
    output: 2VM-sim_quasispecies.fasta
    jobid: 0
    wildcards: sample=2VM-sim

        haploclique --no_singletons --no_prob0 --limit_clique_size=3 --max_cliques=10000 --bam 2VM-sim_aln.bam 2VM-sim_quasispecies
```

Here, we can verify that the `sample` wildcard is replaced by `2VM-sim`. In addition, we see the execution
command for Haploclique.

**Exercise** Run the pipeline on the other provided dataset, namely `5VM-prot_aln.bam`.



### 3.Running V-pipe

The test datasets for the following exercises are located in the `eg_run_vpipe` folder.  Open a terminal,
change into the directory to which you have saved the scripts and test datasets (e.g. inside the VM
`/media/sf_vboxShared/eg_run_vpipe`). 

V-pipe is designed with hierarchically organized data in mind, and it expects the input files to be grouped by
samples. In this context, samples can refer to, e.g. , patient samples or biological replicates of an experiment.
In order to distinguish, different datasets belonging to the same sample, a second level is expected, and it
can refer to, e.g. , sample dates. Below, we show how the working directory for this exercise looks:

```
eg_run_vpipe
├── references
│   ├── HXB2_2253_3869.fasta
│   ├── cohort_consensus.fasta
│   ├── 2VM_msa_2253_3869.fasta
│   └── 5VM_msa_2253_3869.fasta
└── samples
    ├── 2VM-sim
    │   └── 20170904
    │       └── raw_data
    │           ├── 2VirusMix_simulated_R1.fastq.gz
    │           └── 2VirusMix_simulated_R2.fastq.gz
    └── 5VM-prot
        └── 20130626
            └── raw_data
                ├── 5VirusMix_20130626_R1.fastq.gz
                └── 5VirusMix_20130626_R2.fastq.gz
```

V-pipe requires the raw sequencing data in `fastq` format and the reference sequence in `fasta` format. An
initial reference (here, `cohort_consensus.fasta`) is used to obtain a first rough alignment, which will, in
turn, be used for building the profile-HMM (see [ngshmmalign](https://github.com/cbg-ethz/ngshmmalign) 
). If this reference is not provided, an initial
consensus is generated assembling reads *de novo* using software Vicuna[^Vicuna].  A reference sequence (here,
`HXB2_2253_3869.fasta`) is used for pre-filtering reads before running Vicuna, as well as for the alignment
*"lift-over"*[^3]. Additional reference files, namely `2VM_msa_2253_3869.fasta` and `5VM_msa_2253_3869.fasta`,
correspond to the multiple sequencing alignment for true haplotypes. These files are not required *per se*, but
when provided the accuracy of reconstructed haplotypes is assessed.

[^Vicuna]: Xiao Yang, Patrick Charlebois, Sante Gnerre, Matthew G. Coole, Niall J. Lennon, Joshua Z. Levin, James Qu, Elizabeth M. Ryan, Michael C. Zody, and Matthew R. Henn. [De novo assembly of highly diverse viral populations.](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-475) *BMC Genomics*, 13(475), 2012. doi:[10.1186/1471-2164-13-475](https://doi.org/10.1186/1471-2164-13-475)

[^3]: Sequencing reads from different samples are aligned against the profile-HMM. The *"lift-over"* allow the final alignments to be reported with respect to a reference sequence complying to numbering positions conventionally used for downstream analyses. For instance, for HIV-1 the standard reference sequence is HXB2.

**Exercise** Before running V-pipe, let’s have a look at the expected output. We can visually inspect the
order in which rules would be executed using the option `--dag` (see Fig. 2):

``` shell
snakemake -s vpipe.snake -dag | dot -Tpng > dag.png
```

In addition to graphically visualizing the workflow, it is a good practice to check whether the output files
can be created from the inputs, using option `--dryrun`:

``` shell
snakemake -s vpipe.snake --dryrun -p --cores 2
```

By executing a *"dry run"*, one can check if wildcards are used correctly. Below, we show a snippet of the
output displaying the number of jobs to be executed:

```
Job counts:
    count   jobs
    1   all
    2   convert_to_ref
    2   create_simple_initial
    4   extract
    8   gunzip
    2   hmm_align
    2   msa
    2   preprocessing
    4   sam2bam
    27
```

Let’s now proceed and run V-pipe:

``` shell
snakemake -s vpipe.snake -p --cores 2
```

Here, the option `--cores` tells Snakemake the maximum number of cores it can use for rules that support
parallel execution.



## 4. Parameter tuning and customization

The pipeline can be customized through the configuration file `vpipe.config`. The configuration file is a text
file using a basic structure composed of sections, properties and values.

**Exercise** Previously, we ran the pipeline requesting the alignment files as target outputs. In the following,
we will execute rules for reconstructing viral haplotypes. To do so, we will modify the `vpipe.config` file, by
setting the property `global` under the `output` section to `True`:

``` ini
[output]
snv = False
local = False
global = True
```

After modifying the configuration file, try a *"dry run"* to corroborate that rules `haploclique` and `all` are
the only ones to be executed. Below, we show a snippet of the output:

```
Job counts:
    count   jobs
    1   all
    2   haploclique
    3
```

Then, proceed to run V-pipe as before. Newly created files should be located in 
`samples/*/*/variants/global`.

Similarly, the parameters used by the different rules can be tuned. As an exercise, let’s explore some options
of the Haploclique software for global haplotype reconstruction.

**Exercise** Modify the limit on the size of cliques for the haplotype reconstruction. This parameter is used
as a heuristic to reduce run times and memory requirements. Hints: by default, the `clique_size_limit`
option is set to 3. Also, note that you need to remove previously generated output files, otherwise it would
appear as if the target files had already been generated. To do so, you can use rule `haplocliqueclean` as
follows:

``` shell
snakemake -s vpipe.snake -p haplocliqueclean
```

For more details on user-configurable options, please refer to [the wiki](https://github.com/cbg-ethz/V-pipe/wiki/options).



## 5. Result visualization and interpretation

For visualizing the results of the haplotype reconstruction, we have implemented the rule
`haploclique_visualization`. Here, we use multidimensional scaling for visually inspecting the Hamming
distances among reconstructed haplotypes. In addition, for simulations and well-defined lab mixtures, 
reconstructed haplotypes are compared to true haplotypes whenever the corresponding multiple sequence
alignment file is specified.

**Exercise** We will run V-pipe setting rule `haploclique_visualization` as the target rule for each dataset
separately. This is because one of the datasets corresponds to simulated reads from two HIV strains, namely
*"2VM-sim"*, whereas the other dataset was obtained by sequencing a lab-mixture composed of five HIV strains.
Therefore, we must provided different multiple sequence alignment (MSAs) files for the true haplotypes.

Let us start by running V-pipe for analyzing dataset *"5VM-prot"*. This dataset was constructed from a larger
dataset by extracting reads originally mapping to the viral Protease, hence the name. In addition to the
MSA file, we need to specify the starting and ending positions with respect to the reference. Below, we show
a snippet of the configuration file:

``` ini
[haploclique_visualization]
region_start = 0
region_end = 297
msa = references/5VM_msa_2253_3869.fasta
```

In order to execute the pipeline for obtaining the desired plot, open a terminal and type the following:


``` shell
snakemake -s vpipe.snake samples/5VM-prot/20130626/variants/global/quasispecies_plot.pdf -p
```

As output, we should have obtained a tab-separated file containing the mapping of reconstructed haplotypes
to true haplotypes, as well as the multi-dimensional scaling plot.

Now, let’s run the pipeline for dataset *"2VM-sim"*. To do so, you need to modify the configuration file by
providing the proper MSA, as well as the starting and ending position with respect to the reference. This
dataset was generated for a region spanning the viral protease, as well as the reverse transcriptase, which
correspond to a region of length 1617 bases.

We can compare results from this run to results obtained in [section 2.2](#22-wildcards-for-rule-generalization). To do so, change to `eg_haploclique`
folder and type the following:

``` shell
compute_mds -q 2VM-sim_quasispecies -s 0 -e 1917 -r ../eg_run_vpipe/references/2VM_msa_2253_3869.fasta -p quasispecies_plot.pdf -o quasispecies_mapping.tsv
```



## 6. Extending V-pipe with your own tools

V-pipe serves as a modular and extensible backbone, meaning that users can easily customize the workflow
by adding or excluding rules according to their specific requirements. For instance, a different read mapper
or SNV caller can be added for extended functionality.


**Exercise** Write a rule for including the BWA MEM [^bwamem] aligner. Below, we provide a skeleton:

``` python
rule bwa_mem_align:
    input:
        REF = "/path/to/reference",
        R1 = "/path/to/R1.fastq",
        R2 = "/path/to/R2.fastq",
    output:
        SAM = "/path/to/output",
    params:
        BWA = config.applications[’bwa’]
    conda:
        config.bwa_mem_align[’conda’]
    threads:
        config.bwa_mem_align[’threads’]
    shell:
        """
        # 1. index reference
        {params.BWA} index {input.REF}
        # 4. align
        {params.BWA} mem -t {threads} {input.REF} {input.R1} {input.R2} > {output.SAM}
``` 

[^bwamem]: Heng Li. [Aligning sequence reads, clone sequences and assembly contigs with bwa-mem.](https://arxiv.org/abs/1303.3997) ArXiv e-prints,2013

Now, let us assume that you prefer to use BWA as opposed to ngshmmalign. Precedence of rules can be
adjusted using so-called *"local rules"* and setting the rule order. This is exactly how V-pipes decides whether
to use the initial consensus or obtain it *de novo*, before executing rule `hmm_align`.



## 7. Running V-pipe on a LSF cluster

Until now, we have been running V-pipe on our laptops. However, for most applications, running V-pipe
on a local machine may not be efficient. When running V-pipe on a cluster, every rule is executed as an
independent job, but each job is scheduled only after its input files become available. To schedule jobs on a
LSF cluster, you can run V-pipe using option `--cluster`, e.g. :

``` shell
snakemake -s vpipe.snake --cluster "-M {params.mem} -n {threads} -W {params.time} -R \"rusage[mem={params.mem},scratch={params.scratch}]\" -eo {params.lsferrfile} -oo {params.lsfoutfile}" -j 100 --use-conda
```

Using this command, at most 100 jobs will be launched (option `-j`) and resources (e.g memory, time and
threads) will be requested as indicated under the directive `params` for each rule. In addition, by using option
`--use-conda`, Conda environments specific for each rule will be created and all dependencies will be fetched
from the specified channel.



## A. Alternative installation for UNIX-like operating systems

V-pipe integrates various open-source software packages. In order to avoid installation burdens, as well
as overcome incompatibilities due to different software versions, we provide Conda environments for each
module. Conda (see [documentation](https://conda.io/docs/)
) is a cross-platform package manager that is supported by Snakemake (≥3.9.0).

In order to use this feature, we strongly recommend to install Snakemake using Conda. Additionally, we
advise to install Python 3.5 or later versions, as the required Snakemake versions depend on it.


### A.1. Requirements

For installing V-pipe as described below, you will need:

 - (i) conda, and
 - (ii) git
 

### A.2. Installation steps

Before proceeding with the Snakemake installation, check if Conda is installed and if it is in your `PATH`. Type
in a terminal `conda -V`, if conda is installed, you should see something like the following:

    conda 4.2.13

If you need to install Conda, please refer to [the documentation](https://conda.io/docs/install/quick.html)[^4]

[^4]: We recommend to install `miniconda3`, and make sure that the `.bash_profile` file includes the path to `miniconda3`, taking precedence over previous Conda (or Anaconda) installations.

Now, let’s proceed with the installation of Python (≥3.5). If Python is installed (and if it is in your PATH),
you should be able to read the version number, after typing in a terminal: `python --version`. If Python is
not installed or the version is outdated (< 3.5), please install a newer Python version as follows:

``` shell
conda install python=3.6
```

For [installing Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda), type the following:


``` shell
conda install -c bioconda snakemake
```

If Snakemake was succesfully installed, we can have a look at the command line help with `snakemake --help`.
Note that, e.g. if you wish to use the Conda integration, the command-line option `-–use-conda` must be
added to the execution command.

Here, we have used the *bioconda* channel. A Conda channel refers to a repository where Conda looks for
packages. Many bioinformatics tools are available through the bioconda channel and, particularly, most
V-pipe’s software dependencies are included in this repository. Although Conda is a cross-platform tool,
**bioconda packages are not pre-compiled for the native (win32) subsystem Windows**. For this reason, this installation mode
is only recommended for Linux and Mac OS or the [WSL1](https://blogs.msdn.microsoft.com/wsl/2016/04/22/windows-subsystem-for-linux-overview/)
and [WSL2](https://devblogs.microsoft.com/commandline/announcing-wsl-2/) subsystems of Microsoft Windows 10 (a.k.a *Bash for Windows*).

Finally, the pipeline (script, configuration file, conda environments and test data) can be retrieved from
github. In a terminal, change to the directory into which you wish to clone the repository and type:

``` shell
git clone https://github.com/cbg ethz/V pipe.git
```

### Summary

1. Make sure conda is [installed](https://conda.io/docs/install/quick.html)
2. Install Python ≥ 3.5, if needed
3. [Install Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda)
4. Clone [V-pipe github repository](https://github.com/cbg-ethz/V-pipe)

---

Notes and References
