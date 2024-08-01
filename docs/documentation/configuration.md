
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

(configuring-the-workflow)=
## Configuring the workflow

If you have initiated the work directory with `init_project.sh`, you will have a `config.yaml` file in the work directory. This file contains a boilerplate for the configuration for the workflow. All configuration options are described in the schema below.

```{raw} html
:file: config_schema.html
```