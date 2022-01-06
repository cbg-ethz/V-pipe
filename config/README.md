<!-- markdownlint-disable MD013 -->

# Configuring V-pipe

In order to start using V-pipe, you need to provide three things:

 1. Samples in a specific directory structure
 2. _(optional)_ TSV file listing the samples
 3. Configuration file

The utils subdirectory provides [tools](../utils/README.md#samples-mass-importers) that can assist in importing samples files and structuring them.


## Configuration file

The V-pipe workflow is customized using a structured configuration file called `config.yaml`, `config.json` or, for backward compatibility, `vpipe.config` (INI-like format).

This configuration file is a text file written using a basic structure composed of sections, properties and values. When using [YAML](https://yaml.org/spec/1.0/#id2564813) or [JSON](https://www.json.org/json-en.html) format use these languages associative array/dictionaries in two levels for sections and properties. When using the older [INI format](https://docs.python.org/3/library/configparser.html), sections are expected in squared brackets, and properties are followed by corresponding values.

Further more, it is possible to specify additional options on the command line using Snakemake's `--configfile` to pass additional YAML/JSON configuration files, and/or using Snakemake's `--config` to pass sections and properties in a [YAML Flow style](https://yaml.org/spec/1.2.0/#Flow)/JSON syntax.

Here is an **example** of `config.yaml`:

```yaml
general:
  virus_base_config: hiv

input:
  datadir: samples
  samples_file: config/samples.tsv

output:
  datadir: results
  snv: true
  local: true
  global: false
  visualization: true
  QA: true
```

At minimum, a valid configuration **MUST** provide a reference sequence against which to align the short reads from the raw data. This can be done in several ways:

- by using a [_virus base config_](#virus-base-config) that will provide default presets for specific viruses
- by directly passing a reference .fasta file in the section _input_ -> property _reference_ that will override the default

### virus base config

We provide virus-specific base configuration files which contain handy defaults for some viruses.

Currently, the following _virus base config_ are available:

- [hiv](hiv.yaml): provides HXB2 as a reference sequence for HIV, and sets the default aligner to _ngshmmalign_.
- [sars-cov-2](sars-cov-2.yaml): provides NC\_045512.2 as a reference sequence for SARS-CoV-2, sets the default aligner to _bwa_ and sets the variant calling to be done against the reference instead of the cohort's consensus.

### configuration manual

More information about all the available configuration options and an exhaustive list can be found in [config.html](config.html)
or [online](https://htmlpreview.github.io/?https://github.com/cbg-ethz/V-pipe/blob/master/config/config.html).

### legacy V-pipe 1.xx/2.xx users

If you want to re-use your old configuration
from a [legacy V-pipe v1.x/2.x installation](https://github.com/cbg-ethz/V-pipe/wiki/options)
or [sars-cov2 branch](https://cbg-ethz.github.io/V-pipe/tutorial/sars-cov2/#running-v-pipe)
it is possible, if you keep in mind the following caveats:

- The older INI-like syntax is still supported for a `vpipe.config` configuration file.
  - This configuration will be overridden by `config.yaml` or `config.json`,
    you might want to delete those files from your working directory if you are not using them.
- V-pipe starting from version 2.99.1 follows the [Standardized usage](https://snakemake.github.io/snakemake-workflow-catalog/?rules=true) rules of the
  [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=cbg-ethz/V-pipe)
  - This defines a newer [directory structure](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#distribution-and-reproducibility)
    - samples TSV table is now expected to be in `config/samples.tsv`
      (use the section _input_ ->  property _samples_file_ to override).
    - the per sample output isn't written in the same `samples/` directory as the input anymore, but in a separate directory called `results/`
      (use the section _output_ -> property _datadir_ to override).
    - the cohort-wide output isn't written in a different `variants/` directory anymore, but at at the base of the _output datadir_ - i.e by default in `results/`
      (use the section _output_ -> property _cohortdir_ to specify a different path **relative to the output datadir**).
  - Add the following sections and properties to your `vpipe.config` configuration file to **bring back the legacy behaviour**:

```ini
[input]
datadir=samples
samples_file=samples.tsv

[output]
datadir=samples
cohortdir=../variants
```

As of version 2.99.1, only the analysis of viral sequencing data has been
[extensively tested](https://github.com/cbg-ethz/V-pipe/actions/workflows/run_regression_tests.yaml)
and is guaranteed stable.
For other more advanced functionality you might want to wait until a future release.

## samples tsv

File containing sample unique identifiers and dates as tab-separated values.

**Example:** here, we have two samples from patient 1 and one sample from patient 2:

```tsv
patient1    20100113
patient1    20110202
patient2    20081130
```

By default, V-pipe searches for a file named `config/samples.tsv`, if this file does not exist, a list of samples is built by searching the contents of the input datadir.

Optionally, the samples file can contain a third column specifying the read length. This is particularly useful when samples are sequenced using protocols with different read lengths.

## samples

V-pipe expects the input samples to be organized in a two-level directory hierarchy.

- The first level can be, e.g., patient samples or biological replicates of an experiment.
- The second level can be, e.g., different sampling dates or different sequencing runs of the same sample.
- Inside that directory, the sub-directory `raw_data/` holds the sequencing data in FASTQ format (optionally compressed with GZip).

**For example:**

```lang-none
samples
├── patient1
│   ├── 20100113
│   │   └──raw_data
│   │      ├──patient1_20100113_R1.fastq
│   │      └──patient1_20100113_R2.fastq
│   └── 20110202
│       └──raw_data
│          ├──patient1_20100202_R1.fastq
│          └──patient1_20100202_R2.fastq
└── patient2
    └── 20081130
        └──raw_data
           ├──patient2_20081130_R1.fastq.gz
           └──patient2_20081130_R2.fastq.gz
```

The utils subdirectory contain [mass-importers tools](../utils/README.md#samples-mass-importers) to assist you in generating this hierarchy.
