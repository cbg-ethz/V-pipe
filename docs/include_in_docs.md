


### Alternative B: Cov-Spectrum.org

<!-- Most important -->

This is the standard way for us to generate new signatures (outside of new emergent sub-variant that doesn't have enough sequences on Cov-Spectrum yet).

COJAC has quick explanations in its [README](https://github.com/cbg-ethz/cojac/), you can use similar commands to generate full lists of mutations per variants:




```{note}
The above example uses the free and open ENA database. Access to the GISAID database isn't open and requires a token.
```

*[ENA]: European Nucleotide Archive

Then add headers:

<!-- We download this for three variants. What should go at the header of the other variants?  -->
<!-- Where does the source come from? How should I know the metadata of the other variants? -->

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

<!-- INCLUDE this one -->

Emma Hodcroft publishes curated lists of mutation in this directory on Github:
- https://github.com/hodcroftlab/covariants/blob/master/defining_mutations/

As she works using phylogenetic tree, she also flags the reversions. We usually collaborate with her and check each other's mutations lists.

COJAC can then extract a mutation list with:

```bash
curl -O 'https://github.com/hodcroftlab/covariants/raw/master/defining_mutations/23B.Omicron.tsv'
cojac sig-generate --covariants 23B.Omicron.tsv | tee xbb_1_19_mutations_full.yaml
# returns: Error: Missing option '--var' / '--variant'.
```

Finally, add a header to the YAML, with at least a short name, a Pangolineage, a `mut:` section with the mutation list moved there, and a `revert:` section with the reversions.

### Alternative D: UKHSA

<!-- Not important, but mention -->

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


### Interpreting results

Check the amplicons against background on Cov-Spectrum:

```bash
cojac cooc-curate --amplicons results/amplicons.v41.yaml vocs/delta_mutations_full.yaml
```

<!-- This only returns
(not using access keys on API https://lapis.cov-spectrum.org/open/v2) -->

search amplicon which are mostly prevalent in the family searched.

#### e.g.: for XBB*

<!-- This is an example for another strain -->

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