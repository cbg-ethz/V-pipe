# Benchmarking module

Benchmark Quasispecies assembly methods both on the level of SNVs as well as global haplotypes.

## Usage

To run the workflow, execute the following:

```bash
# locally (remove docker part if on linux)
docker run --rm -v $PWD:/foo --workdir=/foo snakemake/snakemake:stable snakemake -prj1 --use-conda

# on cluster
./run_workflow.sh
```

## Adding new methods

To run a new method/tool as part of the benchmark workflow, add a script to `resources/method_definitions/`.
Each script must be classified as either `local` (produces a VCF file) or global (produces a FASTA file) by adding `# GROUP: local` or `# GROUP: global` respectively.
Method dependencies can be specified as comments.
Conda packages can be added by writing `# CONDA: <package name> = <version>`.
Analogously, PIP packages can be added by writing `# PIP: <package name>`.
Multiple packages can be added by repeating these lines.
A conda environment will then be dynamically generated (when running Snakemake with `--use-conda`).


## Configuring your benchmarking study

### Configuration file
`method_list` List of methods that should be executed. Methods should be defined as described above.  
`replicate_count` Number of replicates per line in `params.csv` that should be created.
`haplotype_generation` Either `distance` or `mutation_rate`. In `distance` haplotypes are generated based on distance pattern, in the mode `mutation_rate` haplotypes are based on mutation, deletion and insertion rates indicated in the column `haplos` in `params.csv`.
`params_path` File path to parameter file `params.csv`.
`master_seq_path` Path to fasts file of master sequence that should be used for the generation of the haplotype population.

### Parameter file
This is a `csv`-file with the following columns:
`seq_tech` illumina, nanopore or pacbio
`seq_mode` shotgun, only for `illumina` we also have amplicon, in case of single_amplicon
`seq_mode_param` Parameters for amplicon mode. read_lenght:overlap.
`read_length` Read length
`genome_size` Genome size
`coverage` Coverage  
`haplos`for distance mode: haplos = n_group1,n_group2,d_group12,d_group1,d_group2,freq_dist,freq_param
for mutation mode: haplos = mutation_rate,insertion_rate,deletion_rate,haplotype_pattern
parameters are seperated with "@"
