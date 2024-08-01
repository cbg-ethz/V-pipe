# Running V-pipe

## Local deployment

## Cluster deployment

With snakemake it is possible to submit jobs on a cluster using the batch submission command-line interface of your cluster.

The opensource platform SLURM by SchedMD is one of the popular systems you might find on clusters (Others include LSF, Grid Engine).

The most user friendly way to submit jobs to the cluster is using a special _snakemake profile_.
[smk-simple-slurm](https://github.com/jdblischak/smk-simple-slurm) is a profile that works well in our experience with SLURM (for other platforms see suggestions in [the snakemake-profil documentation](https://github.com/snakemake-profiles/doc)).

```bash
cd vp-analysis/
# download the profile
git clone https://github.com/jdblischak/smk-simple-slurm.git
# edit simple/config.yaml and either comment out the partition and qos or adapt to your local HPC
cat > smk-simple-slurm/simple/config.yaml <<EOT
cluster:
  mkdir -p logs/{rule} &&
  sbatch
    --cpus-per-task={threads}
    --mem={resources.mem_mb}
    --job-name=smk-{rule}-{wildcards}
    --output=logs/{rule}/{rule}-{wildcards}-%j.out
  #--partition={resources.partition}
  #--qos={resources.qos}
default-resources:
  #- partition=<name-of-default-partition>
  #- qos=<name-of-quality-of-service>
  - mem_mb=1000
restart-times: 3
max-jobs-per-second: 10
max-status-checks-per-second: 1
local-cores: 1
latency-wait: 60
jobs: 500
keep-going: True
rerun-incomplete: True
printshellcmds: True
scheduler: greedy
use-conda: True
EOT
cd work_sarscov2/
./vpipe --dry-run --profile ../smk-simple-slurm/simple/ --jobs 100
cd ../..
```

Snakemakes documentation [introduces the key concepts used in profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles).
Check also [the other options for running snakemake on clusters](https://snakemake.readthedocs.io/en/stable/executing/cli.html#CLUSTER) if you need more advanced uses.

### Dependencies downloading on the cluster

In addition, Snakemake has [parameters for conda](https://snakemake.readthedocs.io/en/stable/executing/cli.html#CONDA) that can help management of dependencies:

- using `-conda-create-envs-only` enables to download the dependencies only without running the pipeline itself. This is very useful if the compute nodes of your cluster are not allowed internet access.
- using `--conda-prefix=`_{DIR}_ stores the conda environments of dependencies in a common directory (thus possible to share and re-use between multiple instances of V-pipe).

```bash
cd  vp-analysis/work_sarscov2/
# First download all bioconda dependencies ahead of time
./vpipe --conda-prefix ../snake-envs --cores 1 --conda-create-envs-only
# And then run on the cluster, the compute node will not need to download anything
./vpipe --dry-run --conda-prefix ../snake-envs --profile ../smk-simple-slurm/simple/ --jobs 100
cd ../..
```

When using V-pipe in production environments, plan the installer's `-p` prefix and `-w` working and snakemake's `--conda-prefix` environments directories according to the cluster quotas and time limits.
For example, consider using `${SCRATCH}` and only move the content of the `results/` directory to long-term storage.
