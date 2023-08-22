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
`seq_tech` The following sequencing technologies can be simulated: `illumina`, `nanopore` or `pacbio`.  
`seq_mode` Two modes are possible: `shotgun` is distributing reads randomly along the genome; `amplicon` mode is generating an amplicon scheme based on the parameters given in `seq_mode_param`. Note that the `amplicon` mode is only available with `illumina` as `seq_tech`.  
`seq_mode_param` Parameters for`amplicon` mode. Read length and overlap of reads for the simulation of the amplicon scheme (`read_lenght:overlap` seperated by `:`). It is also possible to set it to `single_amplicon`, in this setting an artifical bed file is generated with one single amplicon across the whole genome. Parameters for amplicon mode.  
`read_length` Read length of generated reads.   
`genome_size` Genome size of the master sequence.    
`coverage` Mean coverage per position in sample.     
`haplos` There are two modes for the haplotype population as specified in the configuration file. Parameters are seperated with "@":  
- For distance mode: haplos = n_group1,n_group2,d_group12,d_group1,d_group2,freq_dist,freq_param   
- or mutation mode: haplos = mutation_rate,insertion_rate,deletion_rate,haplotype_pattern  
