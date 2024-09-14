
# Configuration

(organizing-data)=
## Organizing data

V-pipe expects the input samples to be organized in a **two-level** directory hierarchy.

* The first level can be, e.g., patient samples or biological replicates of an experiment.
* The second level can be, e.g., different sampling dates or different sequencing runs of the same sample.
* Inside that directory, the sub-directory `raw_data` holds the sequencing data in FASTQ format (optionally compressed with GZip). If you use paired end data the files should be named with the suffixes `_R1` and `_R2`.

An example of a directory structure is shown below:

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

By default, V-pipe will look for the samples in the `samples` directory. This can be changed in the configuration file at `input.datadir` see [Configuring the workflow](configuring-the-workflow). 

(setting-up-samplestsv)=
## Setting up `samples.tsv`

Once the samples are organized in the directory structure, the samples need to be listed in a TSV file. This will done automatically after you complete a dry run, e.g.:

```bash
./vpipe --dry-run --cpus 4
```

This will create a `samples.tsv` file in the work directory with the first two columns pre-filled. Of course, you can also create it yourself. In total it can contain four different columns of which the first two are mandatory:

- column 1: first hierarchical level
- column 2: second hierarchical level
- column 3: read length
- column 4: protocol name

In the example above `samples.tsv` would be:

```text
patient1	20100113
patient1	20110202
patient2	20081130
```

Make sure that before you run the pipeline the `samples.tsv` file is correctly filled in. So for example, if you have a read length different then the default 250, either add it as a third column or change the default value in the configuration file at `input.read_length`. 

(specifying-timeline-and-location-information)=
## Specifying timeline and location information

Deconvolution with [Lollipop](https://github.com/cbg-ethz/LolliPop) is based on time-series information. Therefore, you need to provide time information associated with your samples. This can be done in two different ways:

- Providing a `timeline.tsv`
- Providing a `regex.yaml`

```{note}
See [LolliPop's README.md](https://github.com/cbg-ethz/LolliPop#run-the-deconvolution) for more information about configuring the deconvolution.
```

### Provinding a `timeline.tsv`

The file `timeline.tsv` contains the same information as the `samples.tsv` file, but with the addition of the location of the sample. An example for the first few samples of our dataset would be:

```
sample	batch	reads	proto	location_code	date	location
sample1	2021-11-15	251	v41	Ba	2021-11-15	Basel (BS)
sample1	2021-11-16	251	v41	Ba	2021-11-16	Basel (BS)
sample1	2021-11-17	251	v41	Ba	2021-11-17	Basel (BS)
```

```{note}
Note that:

- The timeline tsv contains a header line (`samples.tsv` does not).
- In addition to the first four columns of `samples.tsv`, only `location` and `date` are necessary for LolliPop. The others are optional.
```

### Providing a `regex.yaml`

Often, sample information (including date of sampling) is included in the sample name. We can extract this date information using a regex. V-pipe do this extraction for you, it just needs the regular expression, that you need to provide in a `regex.yaml` file. Here is an example of a regex for the file names we typically use (`PLANT_YEAR_MONTH_DAY_PROPERTIES`):

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

In this case, `regex.yaml` (specified in section  `timeline:`, property `regex_yaml` of the configuration) defines regular expressions that help parse the samples names specified in the above scheme:

```yaml
sample: (?P<location>\d+)_(?P<year>20\d{2})_(?P<month>[01]?\d)_(?P<day>[0-3]?\d)
```

The yaml can contain the following items:
- `sample` and/or `batch`: define regular expressions that are run against the first (and optionally second) column of V-pipe's `samples.tsv`. 
- `datefmt`: [strftime/strptime format string](https://docs.python.org/3/library/datetime.html#strftime-and-strptime-format-codes) to be used on regex named group `date` (e.g.: use `"%Y%m%d"` to parse YYYYMMDD). Specifying `datefmt` is most useful for date formats that don't split nicely into the ` year`, `month`, and `day` regex  named groups: e.g. if your date format uses week number, day of the week, or day of year. In that case, write a regular expression that provides a named-group `date`, and then use, e.g., `%W%w` or `%j` in your ` datefmt`.

The regular expression can contain the following named-groups that are used to build the timeline:

- `location`: this named-group gives the code for the location (e.g.: Ewag's number code in the schema above)
- `year`: year (in `YYYY` or `YY` format. `YY` are automatically expanded to `20YY` --- Yes, I am optimistic with the duration of this pandemic. Or pessimistic with long term use of V-pipe after the turn of century ;-) ).
- `month`: month
- `day`: day
- `date`: an alternative to the year/month/day groups, if dates aren't in a standard format.

Regex are parsed with the [Python regex library](https://pypi.org/project/regex/), and multiple named groups can use the same name. If you expect multiple formats of the file names, you can thus have a construction where you use `|` to give multiple alternative as long as each provide named-groups `location` and either  `year`, `month`, and `day` or `date`:
   
```
(?:(?P<location>\d+)_(?P<year>20\d{2})_(?:(?:(?P<month>[01]?\d)_(?P<day>[0-3]?\d))|(?:R_(?P<repeat>\d+))))|^(?:(?P<location>KLZHCo[vV])(?P<year>\d{2})(?P<month>[01]?\d)(?P<day>[0-3]?\d)(?:_(?P<location_extra>\w+))?)|^(?:(?P<location>B[aA])(?P<BAsam>\d{6})(?:[-_](?P<year>20\d{2})-(?P<month>[01]?\d)-(?P<day>[0-3]?\d))?)
```

(I swear I have personally typed the line above. It has nothing to do with cats walking on my keyboard ฅ^•ﻌ•^ฅ )

### Expanding the location

The short wastewater treatment plant's code (from regex named group `location` in the previous file, or `location_code` in `timeline.tsv`) can be expanded in to the full location name using the file `wastewater_plants.tsv` (this one is specified in the property `locations_table`), e.g.:

```
code    location
10  Zürich (ZH)
16  Genève (GE)
Ba  Basel (BS)
```

(configuring-the-workflow)=
## Configuring the workflow

If you have initiated the work directory with `init_project.sh`, you will have a `config.yaml` file in the work directory. This file contains a boilerplate for the configuration for the workflow. All configuration options are described in the schema below.

```{raw} html
:file: config_schema.html
```