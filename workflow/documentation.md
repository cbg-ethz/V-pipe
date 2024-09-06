# Documentation

As specified in the [usage section](../README.md#usage) of the main README file, these are the steps you need to perform to  use V-pipe.

To configure V-pipe refer to the documentation present in [config/README.md](../config/README.md).

V-pipe expects the input samples to be organized in a [two-level](../config/README.md#samples) directory hierarchy, and the sequencing reads must be provided in a sub-folder named `raw_data`.
Check the utils subdirectory for [mass-importers tools](../utils/README.md#samples-mass-importers) that can assist you in generating this hierarchy.

We provide [virus-specific base configuration files](../config/README.md#virus-base-config) which contain handy defaults for, e.g., HIV and SARS-CoV-2. Set the virus in the general section of the configuration file:
```yaml
general:
  virus_base_config: hiv
```

## Tutorials

If you want to test your new V-pipe installation, we strongly encourage you to check our tutorials which provide real example data.

Tutorials for your first steps with V-pipe for different scenarios are available in the [docs/](../docs/README.md) subdirectory.
