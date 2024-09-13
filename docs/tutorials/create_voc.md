# Prepare Variant of Concern (VOC) data

## COJAC 

The second part of this file (`mut`) can be generated using [COJAC](https://github.com/cbg-ethz/cojac), which queries the Cov-Spectrum database to identify the mutations that are characteristic of each variant. To use COJAC, we need to install it first. We will create a new conda environment called `cowwid-prepare` that contains the necessary tools for preparing the input data. We will use the `mamba` package manager to create the environment and install the required tools.

```bash
# activate the base conda environment
. vp-analysis/*forge*/bin/activate ''

# create the environment
mamba create -n cowwid-prepare -c conda-forge -c bioconda cojac viramp-hub

# deactivate conda
conda deactivate
```

```{note}
If you are using your own conda installation, you can skip `. vp-analysis/*forge*/bin/activate cowwid-prepare`.
```

After installation we use COJAC to generate the yaml files for the delta, omicron BA.1, and omicron BA.2 variants. First activate the `cowwid-prepare` environment:

```bash
# activate the environment 'cowwid-prepare' which contains cojac
. vp-analysis/*forge*/bin/activate cowwid-prepare
```

And create a work directory for our analysis, including the directory where we will store the variant information:

```bash
mkdir -p vp-analysis/work_cowwid/vocs
```

Now we can use COJAC to create the mutation lists for the delta, omicron BA.1, and omicron BA.2 variants:

```bash
cd vp-analysis/work_cowwid/
cojac sig-generate --url https://lapis.cov-spectrum.org/open/v2 --variant B.1.617.2 | tee vocs/delta_mutations_full.yaml
cojac sig-generate --url https://lapis.cov-spectrum.org/open/v2 --variant BA.1 | tee vocs/omicron_ba1_mutations_full.yaml
cojac sig-generate --url https://lapis.cov-spectrum.org/open/v2 --variant BA.2 | tee vocs/omicron_ba2_mutations_full.yaml
```

After creating the yaml files, we need to manually add metadata information
