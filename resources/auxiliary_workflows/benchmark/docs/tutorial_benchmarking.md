---
jupyter:
  jupytext:
    cell_metadata_filter: -all
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.0
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---


# Benchmarking Module Tutorial

The benchmarking module is written using the Snakemake workflow management system. In this tutorial, we explain how to setup your own synthetic benchmark. We will focus on the evaluation of two global haplotype reconstruction methods. However, you can also compare other methods by simply changing the methods executed, and if applicalbe, adapting the performance measures.

The setup consists of two main steps:
1. Specify the data generation (haplotype generation, sequencing read simulation)
2. Choose methods to be included in the study

## Requirements

V-pipe is optimised for Linux or Mac OS systems. Therefore, we recommend users with a Windows system to [install WSL2](https://learn.microsoft.com/en-us/windows/wsl/install) - this is not a full virtual machine but rather a way to run Windows and Linux cooperatively at the same time.

On your system you should have Python2 and Snakemake installed.


## Configuration of your benchmark study

In the directory `V-pipe/resources/auxiliary_workflows/benchmark/config` you find a two files:
- config.yaml  
- params.csv  

See `V-pipe/resources/auxiliary_workflows/benchmark/README.md` for a detailed description.

In the configuration file `config.yaml` the methods that should be executed are specified. Note that the methods in `method_list` should have the same name as the Python file in `resources/method_definitions/` that is executing them. In the example below, we only apply two methods PredictHaplo and CliqueSNV.

```bash
cat <<EOT > ./config/config.yaml
method_list: [predicthaplo, cliquesnv]

replicate_count: 10 # number of replicates generated of for each set of haplotype population parameters

haplotype_generation: mutation_rate # which haplotype population generation method should be used
params_path: config/params.csv  # path to parameter file that is specifying the simulation settings

master_seq_path: ../../../../hiv/HXB2.fasta # path to reference sequence that should be used, if None, a random sequence is produced.

EOT
```

In the parameter file `params.csv`, you can specify the structure of the haplotype population and the sequencing read simulation setup. Each row in the file is a separate sample simulation setup.
An Illumina read sample is simulated consisting of two haplotypes.

```bash
cat <<EOT > ./config/params.csv
seq_tech,seq_mode,seq_mode_param,read_length,genome_size,coverage,haplos
illumina,shotgun,,249,9719,1000,0.1@0.05@0.05@0.6:0.4

EOT
```


## Running your benchmarking study

Before running check what will be executed:
```bash
snakemake --dry-run
```

This command should be executed in the benchmark directory: `V-pipe/resources/auxiliary_workflows/benchmark/`.

To run the full study execute:
```bash
snakemake --use-conda -c1
```
With the command `--use-conda` conda enviroments are installed for each processing step.


## Output directory

The module creates a `results` directory of the following structure:

```text
📁results
├───📁simulated_reads
│   
└───📁method_runs
│   
└───📁haplo_stats
│   
└───📁performance_measures
```

In `simulated_reads` you find the simulated read samples, in `method_runs` are the directories for each sample and method with the corresponding output files of the executed method.
In `haplo_stats` are some summary statistics of the simulated haplotype populations.
In the directory `performance_measures` are the performance comparison tables of the methods.
In our case, where we evaluate the global haplotype reconstruction methods, these are precision, recall and N50 score.
