# Specifying the timeline

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

```
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