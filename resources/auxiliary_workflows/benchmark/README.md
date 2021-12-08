# Quasispecies Reconstruction Benchmark

Benchmark Quasispecies assembly methods both on the level of local as well as global haplotypes.

## Usage

To run the workflow, execute the following:

```bash
# locally
snakemake -prj1 --use-conda

# on cluster
./run_workflow.sh
```

## Adding new methods

To run a new method/tool as part of the benchmark workflow, add a script to `resources/method_definitions/`.
