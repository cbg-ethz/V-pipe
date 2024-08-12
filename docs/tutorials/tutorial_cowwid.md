---
jupyter:
  jupytext:
    cell_metadata_filter: -all
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.5
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

<!-- markdownlint-configure-file { "MD010": { "ignore_code_languages" : [ "tsv", "bash" ] } } -->

# SARS-CoV-2 Wastewater Surveillance Tutorial

This tutorial introduces the type of analysis that we perform as part of [national surveillance program of SARS-CoV-2 variants in wastewater](https://cov-spectrum.ethz.ch/story/wastewater-in-switzerland).
It covers some specific tools that were developped and integrated into V-pipe for the specific analysis of wastewater:

- [COJAC](https://github.com/cbg-ethz/cojac): Integrated into V-pipe, component tool for early detection based on combination of mutations
- [LolliPop](https://github.com/cbg-ethz/LolliPop): Integrated into V-pipe, component tool for kernel-based deconvolution of variants


After reading, you should be able to reproduce the steps necessary to analyse your data
The exact producedure we use in our surveillance program is [documented in the repository COWWID](https://github.com/cbg-ethz/cowwid), and you can refer to that repository for details such as signature files used, and other similar settings.

For the purpose of this Tutorial, we will work with the master branch of V-pipe and use the _sars-cov-2_ virus base config which is adapted for the SARS-CoV-2 virus.

The data we will be using are heavily down-sampled real-world samples that where used in publication:

- Bagutti, Claudia, Monica Alt Hug, Philippe Heim, Laura Maurer Pekerman, Evelyn Ilg Hampe, Philipp Hübner, Simon Fuchs, et al. 2022. "Wastewater Monitoring of SARS-CoV-2 Shows High Correlation with COVID-19 Case Numbers and Allowed Early Detection of the First Confirmed b.1.1.529 Infection in Switzerland: Results of an Observational Surveillance Study." _Swiss Medical Weekly_ 152 (2526): w30202. [doi:10.4414/smw.2022.w30202](https//doi.org/10.4414/smw.2022.w30202)

<!-- Not necessary to mention here, but include elsewhere in the documentation: -->
This set of samples have been sequenced on Illumina sequencers as paired-ends, but V-pipe has been recently adapted to also work on Oxford Nanopore Technologies long reads.


## Requirements

The tutorial assumes that you have installed V-pipe using the [quick install installation documentation](quick-install-v-pipe-and-conda), and that the workflow is setup with the following structure:

```text
vp-analysis
├── V-pipe
├── mambaforge
└── work
```

- `vp-analysis` is the main directory where you have installed V-pipe
- `V-pipe` is the directory with V-pipe's own code
- `mambaforge` has dependencies to start using V-pipe (bioconda, conda-forge, mamba, snakemake)
- `work` is the directory where you have performed the test analysis

### Installing a stand-alone COJAC and viramp-hub

We are going to need to run additional tools for the preparation of data before the analysis.
Let's install them into a separate conda environment:

```bash
# activate the base conda environment
. vp-analysis/*forge*/bin/activate ''
# NOTE: if you use your own conda installation instead, make sure it is active

# create the environment
mamba create -n cowwid-prepare -c conda-forge -c bioconda cojac viramp-hub

# deactivate conda
conda deactivate
```


## Overview of the analysis


### Upstream of analyzing

- getting variants definitions
- setting up parameters related to the mulitplex PCR amplification protocol
- (obviously) getting raw fastq.gz files (or alignments .bam)

### Analysis

- installing and configuring V-pipe
- running COJAC with V-pipe
  - answers: _Which variants are **present**?_
- providing addition information for LolliPop
- running LolliPop with V-pipe
  - answers: _**How much** of the variants are in the mix over time?_

## Signatures for variants

To detect early presence and measure relative abundance of variants, the mutations that occurs on the genome of these variants needs to defined. 
The tools use a custom YAML format to describe variants. The general format of the variant definitions look like this:

```yaml
variant:
  voc: 'VOC-21APR-02'
  who: 'delta'
  short: 'de'
  pangolin: 'B.1.617.2'
  nextstrain: '21A'
source:
- https://github.com/cov-lineages/pango-designation/issues/49
mut:
  …list goes here…
```

```{note}
 - short name in `short` (used for column names internally, etc.)
 - Pangolineage in `pangolin`
 - `mut` contains the list above
 - (the other fields -- Nextstrain, WHO, voc, … -- are optional)
```

There are different possibilities to obtain or produce the necessary YAML files, the exact strategy will depedent on the variant considered (e.g.: widespread vs. new emerging).

### Alternative A: share ours

The variants definitions that we currently use as part of our SARS-CoV-2 variant surveillance in the wastewater are hosted in the repository COWWID in the subdirectory [`voc/`](https://github.com/cbg-ethz/cowwid/tree/master/voc).
Get in touch with us if you have questions, would like to generate new version, etc.
We have access to the GISAID version of Cov-Spectrum.org (consensus sequences are released there earlier, and there can be earlier enough sequences to generate a list for a new variant).


### Alternative B: Cov-Spectrum.org

This is the standard way for us to generate new signatures (outside of new emergent sub-variant that doesn't have enough sequences on Cov-Spectrum yet).

COJAC has quick explanations in its [README](https://github.com/cbg-ethz/cojac/), you can use similar commands to generate full lists of mutations per variants:

<!-- below returned an error -->
  <!-- File "/Users/geertvangeest/Documents/repositories/test-vpipe/vp-analysis/mambaforge/envs/cowwid-prepare/lib/python3.12/json/decoder.py", line 355, in raw_decode
    raise JSONDecodeError("Expecting value", s, err.value) from None
json.decoder.JSONDecodeError: Expecting value: line 1 column 1 (char 0) -->

```bash
# activate the environment 'cowwid-prepare' which contains cojac
. vp-analysis/*forge*/bin/activate cowwid-prepare
# NOTE: if you use your own conda installation, you can simply type `activate cowwid-prepare instead`

cojac sig-generate --url https://lapis.cov-spectrum.org/open/v1 --variant B.1.617.2 | tee delta_mutations_full.yaml
cojac sig-generate --url https://lapis.cov-spectrum.org/open/v1 --variant BA.1 | tee omicron_ba1_mutations_full.yaml
cojac sig-generate --url https://lapis.cov-spectrum.org/open/v1 --variant BA.2 | tee omicron_ba2_mutations_full.yaml
```

```{note}
The above example uses the free and open ENA database. Access to the GISAID database isn't open and requires a token.
```

*[ENA]: European Nucleotide Archive

Then add headers:

```yaml
variant:
  voc: 'VOC-21APR-02'
  who: 'delta'
  short: 'de'
  pangolin: 'B.1.617.2'
  nextstrain: '21A'
source:
- https://github.com/cov-lineages/pango-designation/issues/49
mut:
  …list goes here…
```


### Alternative C: Covariants.org

Emma Hodcroft publishes curated lists of mutation in this directory on Github:
- https://github.com/hodcroftlab/covariants/blob/master/defining_mutations/

As she works using phylogenetic tree, she also flags the reversions. We usually collaborate with her and check each other's mutations lists.

COJAC can then extract a mutation list with:

```bash
curl -O 'https://github.com/hodcroftlab/covariants/raw/master/defining_mutations/23B.Omicron.tsv'
cojac sig-generate --covariants 23B.Omicron.tsv | tee xbb_1_19_mutations_full.yaml
```

Finally, add a header to the YAML, with at least a short name, a Pangolineage, a `mut:` section with the mutation list moved there, and a `revert:` section with the reversions.

### Alternative D: UKHSA

UK HSA publishes their own variant definitions in this repo:
 - https://github.com/ukhsa-collaboration/variant_definitions

*[HSA]:  Health Security Agency
*[UKHSA]:  United Kingdom's Health Security Agency

```{note}
These definitions are geared toward the typing of consensus sequences and aren't exhaustive. In our experience, due to the dispersion nature of wastewater sequencing, exhaustive list usually perform better, as a smaller curated subset like UKHSA's might all fall on a drop outs.
```

It's possible to convert their YAML format into COJAC's by using:

```bash
phe2cojac --shortname 'om2' --yaml voc/omicron_ba2_mutations.yaml variant_definitions/variant_yaml/imagines-viewable.ym
```

```{note}
a short name needs to be passed on the command line, the rest of the header is generated out of information available in the converted YAML.
```

### Checking background levels of mutations

Once you have a YAML, check the mutation list with `cojac cooc-curate`

## Installation

<!-- This requires explanation -->
We need to request a specific branch:

```bash
curl -O 'https://raw.githubusercontent.com/cbg-ethz/V-pipe/ninjaturtles/utils/quick_install.sh'
bash quick_install.sh -b ninjaturtles -p vp-analysis -w work
```

- use `-b` option for branch

will produce:

```text
📁vp-analysis
├───📁V-pipe      # V-pipe checked out from Github
├───📁mambaforge  # bioconda + conda-forge + mamba + Snakemake
└───📁work        # work directory
```

<!-- What's this? -->

curl 'https://polybox.ethz.ch/index.php/s/ZumerimckhgCRYB/download' -o cowwid-tutorial.tar.xz

## Getting data into V-pipe

<!-- Where do we get this data? -->
What we need is a `samples.tsv` and the specific 2-level hierarchy that V-pipe uses (1-level should also be okay):

```text
📁samples
├───📁sample1
│   └───📁sequencingdate1
│       └───📁raw_data
│           └───🧬reads.fastq
└───📁sample2
    ├───📁sequencingdate1
    │   └───📁raw_data
    │       └───🧬reads.fastq
    └───📁sequencingdate2
        └───📁raw_data
            └───🧬reads.fastq
```

### Alternative A: Getting .fastq.gz files

There are also tools useful for automatically importing files:

<!-- What does this do? -->
<!-- Also not documented in the script -->

```bash
../V-pipe/utils/sort_samples_dumb -h
../V-pipe/utils/sort_samples_dumb -b 20230331_HN3YHDRX2 -f ww_benchmark/samples/ -t samples.imported.tsv -o samples/
```

Important for now, add column 4 with proto:

```bash
gawk 'BEGIN{OFS="\t"};{print $0, "v41"}' samples.imported.tsv |tee samples.tsv
```

> That's how the automation of the current clinical+wastewater surveillance gets the data into V-pipe

### Alternative B: Getting already aligned .bam files

e.g. importing output generated by Artic's own workflow. (No tool available yet)

```text
results
├──📁sample1
│  ├──📁20100113
│  │  └──📁alignments
│  │     └──REF_aln.bam
│  └──📁20110202
│     └──📁alignments
│        └──REF_aln.bam
└─📁sample2
  …etc…
```

- if they are already trimmed, name them `REF_aln_trim.bam` instead.

> That's how currently we pass data for the wastewater surveillance between the production V-pipe (which runs most of the surveillance) and the special branch that runs COJAC and LolliPop.

## Configuration

sources:
 - the quick introduction in the file [config/README.md](https://github.com/cbg-ethz/V-pipe/tree/ninjaturtles/config)
 - the full manual in config/config.html (doesn't work online, only locally from your disk).

```yaml
general:
    virus_base_config: 'sars-cov-2'
    primers_trimmer: samtools
    # for Oxford nanopore
    aligner: minimap
    reprocessor: skip

input:
    datadir: samples/
    samples_file: samples.tsv
    # for Oxford nanopore
    paired: false
    # generated with COJAC (or obtained from us)
    variants_def_directory: references/voc/

output:
    datadir: results/

    trim_primers: true
    snv: false
    local: false
    global: false
    visualization: false
    diversity: false
    QA: false
    upload: false
    dehumanized_raw_reads: false
    # note no wastewater output flag for now, rules called explicitly

# for Oxford nanopore
minimap_align:
    preset: 'map-ont'

# if dates and location are extracted from sample names:
timeline:
    regex_yaml: regex.yaml
    locations_table: wastewater_plants.tsv

deconvolution:
    threads: 8
    # this file corresponds to the parameters used now on our curves:
    # (provided by us)
    deconvolution_config: deconv_bootstrap_cowwid.yaml
    # file that specifies which variant are present at which time point, as determined by looking at COJAC's results
    # done manually by user
    variants_dates: var_dates.yaml
    # automatically generated
    variants_config: results/variants_pangolin.yaml
```

- The `deconvolution_config` parameters points to presets for the algorithm generating the curves.
  - LolliPop has ready-to-use YAMLs directory [`presets/`](https://github.com/cbg-ethz/LolliPop/tree/main/presets)
  - Currently we use, [`deconv_bootstrap_cowwid.yaml`](https://github.com/cbg-ethz/LolliPop/tree/main/presets/deconv_bootstrap_cowwid.yaml) in production.
  - We will eventually switch to [`deconv_linear_logit_quasi_strat.yaml`](https://github.com/cbg-ethz/LolliPop/blob/main/presets/deconv_linear_logit_quasi_strat.yaml) (presented in the [LolliPop pre-print](https://www.medrxiv.org/content/10.1101/2022.11.02.22281825v1)).
    It's much faster as it doesn't rely on bootstrapping, bug currently confidence intervals have some instability at low concentrations.
- The above example also demonstrates a few options for running on Oxford Nanopore data (Setting the aligned to _minimap_ instead of SARS-CoV-2's default _bwa_, single reads, etc.)

## Run COJAC

COJAC helps answer the question:
- is a given variant **present** in the water?

### Commands

This command will run all the COJAC processing, and generate a report-like CSV for protocol Artic V4.1:

```bash
./vpipe --cores 8 allCooc results/cohort_cooc_report.v41.csv
```

(And this command will generate just the amplicons list for Artic V4.1, to checks them before running the rest:)

```bash
./vpipe --cores 8 results/amplicons.v41.yaml
```

### Interpreting results

Check the amplicons against background on Cov-Spectrum:

```bash
cojac cooc-curate --amplicons results/amplicons.v41.yaml
```

search amplicon which are mostly prevalent in the family searched.

#### e.g.: for XBB*

```bash
cojac cooc-curate --amplicons amplicons.v41.yaml references/voc/xbb_mutations_full.yaml | tee xbb_amplicon_curate.ansi
```

Despite **19326G** being exclusive to XBB, it's not close to any other mutation so it's not possible to look for it as a combination of multiple cooccurrences (otherwise, see ["Other situations" below](other situations)), **BUT** that part is interesting:

> 75_omxbb[22577CA,22599C,22664A,22674T,22679C,22686T,22688G,22775A]: ***XBB*=0.85**, BJ.1=0.49, BA.2.10.1=0.00
> ***76_omxbb[22775A,22786C,22813T,22882G,22895CC,22898A,22942G,22992A,22995A,23013C,23019C]: *XBB*=0.85**, BJ.1=0.01
> 77_omxbb[22992A,22995A,23013C,23019C,23031C,23055G,23063T,23075C]: ***XBB*=0.94**, BM.1.1.1=0.75, BM.4.1=0.02

**Note** The **_emphasis_** is on the variant family considered.

On the ARTIC v4.1 amplicon number 76, despite none of the mutations being exclusive to XBB, this peculiar _combination_ is exclusive to XBB according to Cov-Spectrum.
Thanks to 22664A and 22895CC being somewhat more frequent in XBB (but also BJ.1), and the other mutation being most frequent in _different_ variants (e.g.: 22942G and 23019C are _never found_ in BJ.1)

Then one can look at content of `results/cohort_cooc_report.v41.tsv`, or use `cojac cooc-colormut` with `results/cohort_cooc.v41.yaml`.

> **Tip** You can also edit a subset of amplicon.yaml and use that when running cojac display tools.

### Other situations

Sometimes there are no clear amplicons for detecting a variant.
Other strategies including tracking multiple single mutations.

The option `mincooc` in the section `amplicon:` of the configuration controls how many mutation cooccurrences at minimum are considered per amplicon. By default it is 2 (consider amplicons carrying a duplet of mutations), but by lowering to 1 it is also possible to search for singleton mutations. Remember that single mutations aren't very informative, so try to combine information from several to increase confidence.

(TO BE DOCUMENTED LATER)

## LolliPop

LolliPop helps answer the question:
 - In which **relative proportions** are the variants in the water?

### Variants per dates

The file `var_dates.yaml` (specified in the section `deconvolution:` of the configuration) gives information to LolliPop which variants to run deconvolution on which time period.
We use it to inform LolliPop which variants we know to be present:
- based on COJAC output
- based on other sources (e.g.: clinical case detection)
- based on general information (NextStrain's _Molecular clock_, the date of discovery of a new variant)

```yaml
var_dates:
  '2022-08-15':
  - BA.4
  - BA.5
  - BA.2.75
  #- BA.2.75.2
  - BQ.1.1
  '2022-11-01':
  - BA.4
  - BA.5
  - BA.2.75
  - BQ.1.1
  - XBB
```

### Timeline

Because it doesn't treat each sample separately, but considers them as a time series and leverages this to compensate for dispersion (it relies on a _kernel-based_ deconvolution), LolliPop needs additional information (locations and sampling dates) which is not provided in the standard V-pipe `samples.tsv` yet.

The timeline file is a file in TSV format and adds extra columns:

```tsv
sample	batch	reads	proto	location_code	date	location
A1_05_2023_04_12	20230428_HNG5MDRX2	250	v41	5	2023-04-12	Lugano (TI)
A2_10_2023_04_13	20230428_HNG5MDRX2	250	v41	10	2023-04-13	Zürich (ZH)
A3_16_2023_04_14	20230428_HNG5MDRX2	250	v41	16	2023-04-14	Genève (GE)
…
```

- The extra columns _location_ and _data_ are the **necessary** one.
- Columns _sample_, _batch_, _reads_ and  _proto_ are simply the fist four columns of `samples.tsv`
  - V-pipe only uses column _sample_ and _batch_ for now.

#### Alternative A: regex.yaml

The way we do it is that we have a fixed naming scheme:

```text
┌──────────────── Wastewater Treatment Plant:
│                  05 - CDA Lugano
│                  10 - ARA Werdhölzli in Zurich
│                  12 - STEP Vidy in Lausanne
│                  17 - ARA Chur
│                  19 - ARA Altenrhein
│                  25 - ARA Sensetal
│  ┌───────────── Date
│  │          ┌── Sample properties
┴─ ┴───────── ┴─
09_2020_03_24_B
10_2020_03_03_B
10_2020_03_24_A
10_2020_04_26_30kd
```

so the file `regex.yaml` (specified in section  `timeline:`, property `regex_yaml` of the configuration) defines regular expressions that help parse the samples names specified in the above scheme:

```yaml
sample: (?P<location>\d+)_(?P<year>20\d{2})_(?P<month>[01]?\d)_(?P<day>[0-3]?\d)
```

- `sample` (and optionally `batch`) define regular expressions that are run against the first (and optionally second) column of V-pipe's `samples.tsv`. They define the following named-groups
  - `location`: this named-group gives the code for the location (e.g.: Ewag's number code in the schema above)
  - `year`: year (in `YYYY` or `YY` format. `YY` are automatically expanded to `20YY` --- Yes, I am optimistic with the duration of this pandemic. Or pessimistic with long term use of V-pipe after the turn of century ;-) ).
  - `month`: month
  - `day`: day
  - `date`: an alternative to the year/month/day groups, if dates aren't in a standard format.
  - regex are parsed with the [Python regex library](https://pypi.org/project/regex/), and multiple named groups can use the same name.
    You can thus have a construction where you use `|` to give multiple alternative as long as each provide named-groups `location` and either  `year`, `month`, and `day` or `date`:
    ```regex
    (?:(?P<location>\d+)_(?P<year>20\d{2})_(?:(?:(?P<month>[01]?\d)_(?P<day>[0-3]?\d))|(?:R_(?P<repeat>\d+))))|^(?:(?P<location>KLZHCo[vV])(?P<year>\d{2})(?P<month>[01]?\d)(?P<day>[0-3]?\d)(?:_(?P<location_extra>\w+))?)|^(?:(?P<location>B[aA])(?P<BAsam>\d{6})(?:[-_](?P<year>20\d{2})-(?P<month>[01]?\d)-(?P<day>[0-3]?\d))?)
    ```
    (I swear I have personally typed the line above. It has nothing to do with cats walking on my keyboard ฅ^•ﻌ•^ฅ ).
- `datefmt`: [strftime/strptime format string](https://docs.python.org/3/library/datetime.html#strftime-and-strptime-format-codes) to be used on regex named group `date` (e.g.: use `"%Y%m%d"` to parse YYYYMMDD).
  - This is most useful for date formats that don't split nicely into the ` year`, `month`, and `day` regex  named groups: e.g. if your date format uses week number, day of the week, or day of year.
    In that case, write a regular expression that provides a named-group `date`, and then use, e.g., `%W%w` or `%j` in your ` datefmt`.

The short wastewater treatment plant's code (from regex named group `location` in the previous file) is then expanded in to the full location name using the file `wastewater_plants.tsv` (this one is specified in the property `locations_table`), e.g.:

```tsv
code    location
10  Zürich (ZH)
16  Genève (GE)
Ba  Basel (BS)
```

You need to adapt this procedure to your needs.
Do not hesitate to contact us and to check the Timeline section of the exhaustive configuration manual in your (locally on your hard-drive: config/config.html).
It is also possible to use other schemes (e.g.: sequencing batch _is_ the sampling date, using dates in different format -- e.g. week number -- etc.)

#### Alternative B: providing your own timeline.tsv

It is also possible to write and provide your own file.

This can either be done prior to starting V-pipe -- e.g. an external software could query your LIMS' database and add the necessary column to sample.tsv in order to generate the table -- specify the location of this output in section `tallymut:` property `timeline_file`.

Or by heavily customizing the timeline rule -- e.g. using the `timeline:` section, property `script` to run your own script instead of V-pipe's official regex-based extractor.

### Others

There two last files controlling LolliPop, but those usually won't require much attention:
- `deconvolution_config` points to presets describing the algorithm generating the curves -- we simply use the [`presets/`](https://github.com/cbg-ethz/LolliPop/tree/main/presets) available from LolliPop repository.
- `variants_config` gives additional information about how to process the variants.
  - at minimum, it should contain a section `variants_pangolin:` mapping _short names_ used in various files back to the full Pangolineages used in the results. V-pipe will automatically generate one (`results/variants_pangolin.yaml`) and use it.
  - otherwise, other sections can be added to specify only a subset of locations, start date, end date, etc.

See [LolliPop's README.md](https://github.com/cbg-ethz/LolliPop#run-the-deconvolution) for more information about configuring the deconvolution.

### Run LolliPop

This command will run the deconvolution:

```bash
./vpipe --cores 8 deconvolution
```

### Interpret results

V-pipe will output two files using LolliPop:
- `results/deconvoluted.tsv.zst`: a Zstandard-compressed table that can be directly reads into pandas for further processing.
- `results/deconvoluted_upload.json`: curves in a JSON format that can be used to upload to dashboard.

For further examples,  in repository COWWID, see files:
- [ww_cov_uploader_V-pipe.ipynb](https://github.com/cbg-ethz/cowwid/tree/master/ww_cov_uploader_V-pipe.ipynb) -- will display the curves in the Notebook, before uploading them to [Cov-Spectrum.org](https://cov-spectrum.ethz.ch/story/wastewater-in-switzerland) and [BAG/FOPH](https://www.covid19.admin.ch/en/epidemiologic/waste-water).
- [DeconvolutionPrediagnostics.ipynb](https://github.com/cbg-ethz/cowwid/tree/master/DeconvolutionPrediagnostics.ipynb) -- help diagnose problems related to data (drop-out affecting mutations) and/or signatures (too much similarity) that could affect the deconvolution.
