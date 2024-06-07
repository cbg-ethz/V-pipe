<!-- markdownlint-disable MD013 -->
<!-- markdownlint-configure-file { "MD010": { "ignore_code_languages" : [ "tsv" ] } } -->
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
  In addition, a look-up for the recent versions of ARTIC protocol is provided; this makes it possible to set per-sample protocol in the sample table, and to turn on amplicon trimming (see [amplicon protocols](#amplicon-protocols)).
- [rsvb](rsvb.yaml) config file for Human respiratory syncytial virus B (RSV-B), used to process Illumina RSV samples.
- [h3n2_ha](h3n2_ha.yaml) config-file used for the analysis of H3N2 segment HA from wastewater Illumina data on SRA available through the SRA Run accession: [SRP385331](https://www.ebi.ac.uk/ena/browser/view/PRJNA856656)
- [drosophila_c_virus](drosophila_c_virus.yaml) configuration used for the analysis of drosphila C virus (DCV) Illumina samples in Lezcano et al., Virus Evolution, 2023, doi:[10.1093/ve/vead074](https://doi.org/10.1093/ve/vead074), NCBI BioProject accession number [PRJNA993483](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA993483)
- [](herpes_simplex_virus_2.yaml) config-file used for the analysis of herpes simplex virus 2 (HSV-2) Illumina samples in Lezcano et al., Virus Evolution, 2023, doi:[10.1093/ve/vead074](https://doi.org/10.1093/ve/vead074). Analysed sample is from  López-Muñoz AD, Rastrojo A, Kropp KA, Viejo-Borbolla A, Alcamí A. "Combination of long- and short-read sequencing fully resolves complex repeats of herpes simplex virus 2 strain MS complete genome". _Microb Genom._ 2021 Jun;7(6). Sample accession number [ERR3278849](https://www.ebi.ac.uk/ena/browser/view/ERR3278849)
- [polio](polio.yaml) config-file used for the analysis of poliovirus MinION samples. sample accession number: [ERR4027774]](https://www.ebi.ac.uk/ena/browser/view/ERR4027774) (Shaw et al., 2020, DOI: https://doi.org/10.1128/jcm.00920-20)
- [mpxv](mpxv.yaml) Monkey pox virus

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
patient1	20100113
patient1	20110202
patient2	20081130
```

By default, V-pipe searches for a file named `config/samples.tsv`, if this file does not exist, a list of samples is built by searching the contents of the input datadir.

### read-lenght

The samples' read-length is used for critical steps of the pipeline (e.g.: quality filtering). Different possibilities are available to set its value:

- by default, V-pipe expects a read-length of 250bp
- this default can be globally overridden in the configuration file in section _input_ -> property _read_length_
  ```yaml
  input:
    read_length: 150
  ```
- the samples TSV file can contain an optional third column specifying the read length.
  This is particularly useful when samples are sequenced using protocols with different read lengths.
  ```tsv
  patient1	20100113	150
  patient1	20110202	200
  patient2	20081130	150
  ```

  The utils subdirectory contain [mass-importers tools](../utils/README.md#samples-mass-importers) that can generate this third column while importing samples.

### amplicon protocols

Samples can be the result of PCR amplification. This can require some additional processing, e.g., primers might need trimming:

```yaml
output:
  trim_primers: true
```

In order to complete these steps, additional information needs to be provided, e.g., a BED file describing the primers to be trimmed.

- This can be specified globally with several properties in the configuration file in section _input_:
  ```yaml
  input:
    primers_bedfile: references/primers/SARS-CoV-2.primer.bed
    inserts_bedfile: references/primers/SARS-CoV-2.insert.bed
  ```

- The samples TSV file can contain an optional fourth column specifying the protocol:
  - When different samples have been processed with different library protocols, a lookup table with per-protocol specific (primers bed and fasta), can be provided in a YAML file.
    `references/primers.yaml`:
    ```yaml
    v41:
      name: SARS-CoV-2 ARTIC V4.1
      inserts_bedfile: references/primers/v41/SARS-CoV-2.insert.bed
      primers_bedfile: references/primers/v41/SARS-CoV-2.primer.bed
    v4:
      name: SARS-CoV-2 ARTIC V4
      inserts_bedfile: references/primers/v4/SARS-CoV-2.insert.bed
      primers_bedfile: references/primers/v4/SARS-CoV-2.primer.bed
    v3:
      name: SARS-CoV-2 ARTIC V3
      inserts_bedfile: references/primers/v3/nCoV-2019.insert.bed
      primers_bedfile: references/primers/v3/nCoV-2019.primer.bed
    ```  
  - in the configuration file, this look-up can be then specified in section _input_ option _protocols_file_:
   `config/config.yaml`:
    ```yaml
    input:
      protocols_file: references/primers.yaml
    ```
  - The short name can now be referenced in the fourth column samples TSV table file:
    `config/samples.tsv`:
    ```tsv
    sample_a	20211108	250	v3
    sample_b	20220214	250	v4
    ```

  This is useful if multiple different amplicon schemes have been used of the lifetime of a long-running project, as new variants appear over time with SNVs that require adapting amplicons.

- [_virus base config_](#virus-base-config) can provide some defaults for either above
  e.g.: sars-cov-2 provides BED files for ARTIC v3, v4 and v4.1


## samples

V-pipe expects the input samples to be organized in a two-level directory hierarchy.

- The first level can be, e.g., patient samples or biological replicates of an experiment.
- The second level can be, e.g., different sampling dates or different sequencing runs of the same sample.
- Inside that directory, the sub-directory `raw_data/` holds the sequencing data in FASTQ format (optionally compressed with GZip).

**For example:**

```text
📁samples
├──📁patient1
│  ├──📁20100113
│  │  └──📁raw_data
│  │     ├──🧬patient1_20100113_R1.fastq
│  │     └──🧬patient1_20100113_R2.fastq
│  └──📁20110202
│     └──📁raw_data
│        ├──🧬patient1_20100202_R1.fastq
│        └──🧬patient1_20100202_R2.fastq
└──📁patient2
   └──📁20081130
      └──📁raw_data
         ├──🧬patient2_20081130_R1.fastq.gz
         └──🧬patient2_20081130_R2.fastq.gz
```

The utils subdirectory contain [mass-importers tools](../utils/README.md#samples-mass-importers) to assist you in generating this hierarchy.
