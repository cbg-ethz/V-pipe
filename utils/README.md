<!-- markdownlint-disable MD013 MD010 -->

# Utils

This directory contains several command-line utilities that can help in some ancillary tasks to get V-pipe running.

## Quick installer

`quick_install.sh` is a script that can assist deploying V-pipe: It can automatically download and install bioconda, snakemake and fetch V-pipe from the repository.

It is possible to directly run it from the web with the following commands:

```bash
curl -O 'https://raw.githubusercontent.com/cbg-ethz/V-pipe/master/utils/quick_install.sh'
bash quick_install.sh -w working
cd ./working/
```

- This will download and install bioconda, snakemake and V-pipe in the current directory (use option `-p` to specify another directory).
- The installed version by default is the `master` git branch
  - It is possible to specify another git branch or a tag using the `-b` options, e.g. `-b v2.99.1`.
  - Alternatively, the `-r` option will download and install the `.tar.gz` tar-ball package of a specific version.
  - See the [release page](https://github.com/cbg-ethz/V-pipe/releases) for available tags and releases.
- using `-w` will create a working directory and populate it.
  - It will copy over a default `config.yaml` that you can [edit to your liking](../config/README.md)
  - And it will create a handy `vpipe` short-cut script to invoke `snakemake`:
```bash
cd working
# edit config.yaml and provide samples/ directory
./vpipe --jobs 4 --printshellcmds --dry-run
```

> **Tips**:
> To create and populate other new working directories, you can call `init_project.sh` from within the new directory:
>
> ```bash
> mkdir -p working_2
> cd working_2
> ../V-pipe/init_project.sh
> ```

Available options:

```console
# ./quick_install.sh -h
usage: ./quick_install.sh [options]
options:
-f           force overwriting directories
-p PREFIX    prefix directory under which to install V-pipe
             [default: current directory]
-b BRANCH    install specified branch of V-pipe
             [default: master]
-r RELEASE   install specified release package of V-pipe
-w WORKDIR   create and populate working directory
-m           only minimal working directory
-h           print this help message and exit"
```


## Samples mass-importers

Here are also some tools that can help mass-importing samples into V-pipe.
They will search for `.fastq.gz` files, put them in the [two-level hierarchy that V-pipe expects in  `samples/`](../config/README.md#samples) and generate the corresponding `samples.tsv`. By default, hard-links will be used to save space and avoid duplicating these large sequencing files.

Whichever of these tools is most suitable for you depends on how you receive the files from the sequencing lab.

- [`sort_samples_dumb`](#sort_samples_dumb) is for loose collection of FASTQ files
- [`sort_samples_demultiplexstats`](#sort_samples_demultiplexstats) and [`sort_samples_jobinfo`](#sort_samples_jobinfo) are better suited when additional information has been provided by the base calling and demultiplexing software.
  - these two also support _patch maps_ TSV files, helping to rename the samples from the name used in the original sample sheet in the LIMS (laboratory's information management system) to something more flexible.
    For example, to remap simple sequential numbers to longer names, use the following patch map TSV:

```tsv
1	05_2021_02_01
2	05_2021_02_02
3	05_2021_02_03
4	05_2021_02_04
5	05_2021_02_05
268	05_2021_02_06
269	05_2021_02_07
270	05_2021_02_08
271	05_2021_02_09
272	05_2021_02_10
273	05_2021_02_11
274	05_2021_02_12
275	05_2021_02_13
```

### sort_samples_dumb

`sort_samples_dumb` is useful when labs are providing the sequencing data simply as a loose collection of FASTQ files.

Example of usage:

```console
V-pipe/utils/sort_samples_dumb -f download/ -t working/samples.tsv -o working/samples -b 20210110
```

Where:

- `-f` : specifies the main directory containing the downloaded the `.fastq.gz` files.
  - all its subdirectories will be searched recursively
- `-t` : specifies the `samples.tsv` file to create
- `-o` : specifies the output directory where to store the files
- `-b` : is the directory to use for the second level (see V-pipe tutorials), e.g.: dates
  - the first level is usually the patient or sample name and `sort_samples_dumb` will attempt to guess it from the file names.
  - the second level is usually the sampling date or sequencing batch, but `sort_samples_dumb` has no simple way to guess that.

Running this command will immediately copy over the `.fastq.gz` files (using hard links) in the working directory and generate the `samples.tsv`. After checking the content, you should be able to run V-pipe.

Available options:

```console
./sort_samples_dumb -h
Usage: ./sort_samples_dumb -f <DIR> -b <BATCH> [-l <LEN>] [-L {''|--link|--symbolic-link|--reflink}] -o <DIR> -m <MODE>
	-f : directory containing .fastq.gz files
	-b : batch name to use for 2nd level (e.g.: date)
	-l : read lenght (default: autodetect)
	-L : link parameter to pass to cp when copying (default: --link)
	-t : tsv file (default: samples.<BATCH>.tsv)
	-T : do not truncate (empty) the file before starting
	-g : store list in .tsv.staging instead and only rename into final .tsv if successful
	-D : sample have duplicates (e.g.: across lanes)
	-p : prefix to prepend to fastq files (e.g.: for fusing runs)
	-s : suffix to append to fastq files (e.g.: for fusing runs)
	-o : output directory
	-m : POSIX mode parameter to pass to mkdir (e.g.: 0770)
```

Of special interest:

- `-l`: sets the read length instead of trying to autodetect it with an _awk_ script.
- `-D`: when multiple files have the same name, the importer will group them for merging as the same sample (e.g.: lane-duplicates).


### sort_samples_demultiplexstats

`sort_samples_demultiplexstats` can help if the lab provides the DemultiPlex Stats files generated by Illumina's `bcl2fastq` version 2.x software. This tool will use the information provided in the JSON file to match files to the respective samples they came from.

The general syntax is:

```console
V-pipe/utils/sort_samples_demultiplexstats --statsdir downloads/Demultiplex --fastqdir downloads/RawReads --qcdir downloads/FastQC --outdir working/samples
```

Where:

- `--statsdir` : is the directory containing the file `Stats/Stats.json` produced by `bcl2fastq`.
- `--outdir` : directory where to create the output
optionally:
- `--fastqdir`: is the directory containing the `.fastq.gz` files if they are not in the same directory as stats.
- `--qcdir`: is the directory with the FastQC's quality checks if provided by the lab.
- `--noempty`: if the lab deleted the empty (0 reads) `.fastq.gz` files.

This command is the **first step** (analysing `Stats.json`):

- It will generate a samples TSV file in the output directory: `working/samples/samples.`_{some date}_`.tsv`.
  You can copy the content of this file into your `samples.tsv`
- It will generate a file `working/samples/movedata.sh`, you can run it with:
  ```console
  bash working/samples/movedata.sh
  ```
  and that file will in turn perform the **second step**: hard-linking all the files from `download/RawReads` into `working/samples/`_{sample name}_`/`_{some date}_.
- once the second step has been performed, you should be able to run V-pipe.

Available options:

```console
./sort_samples_demultiplexstats  -h
usage: sort_samples_demultiplexstats [-h] -S DIR [-f DIR] [-q DIR] [-o DIR] [-m MODE] [-L CPLINK] [-s] [-a] [-n] [-p TSV]

Uses bcl2fastq's demultiplexing stats as metadata to organise samples

optional arguments:
  -h, --help            show this help message and exit
  -S DIR, --statsdir DIR
                        directory containing 'Stats/Stats.json'
  -f DIR, --fastqdir DIR
                        directory containing .fastq.gz files if different from above
  -q DIR, --qcdir DIR   if set, import FastQC's _fastqc.html files from there
  -o DIR, --outdir DIR  output directory
  -m MODE, --mode MODE  POSIX file access mode to be passed to mkdir
  -L CPLINK, --linking CPLINK
                        parameter to pass to `cp` for linking files instead of copying their data
  --force               Force overwriting any existing file when moving
  -s, --summary         Only display a summary of datasets, not an exhaustive list of all samples
  -a, --append          Append to the end of movedatafiles.sh, instead of overwritting (use when calling from an external combiner wrapper)
  -g, --staging         Write samples list in .tsv.staging and only rename them to the final .tsv at the end of movedatafiles.sh if there were no errors.
  -n, --noempty         skip fastq.gz files with bad yield (0 reads)
  -p TSV, --patchmap TSV
                        patchmap file to rename samples
```


### sort_samples_jobinfo

`sort_samples_jobinfo` is similar in concept to [`sort_samples_demultiplexstats`](#sort_samples_demultiplexstats), but it relies on the `CompletedJobInfo.xml` and `SampleSheetUsed.csv` files generated by the Illumina Analysis Software on Windows, if those are files are provided by the lab, and will try to match samples listed therein to  `.fastq.gz` files.

The general syntax is:

```console
V-pipe/utils/sort_samples_jobinfo --sourcedir=downloads/20210528_061936 --outdir=working/samples
```

Where:

- `--sourcedir` : is the directory, usually named _{yymmdd}_`_`_{hhmmss}_ and found inside the subdirectory `Alignment_1` of the Illumina's run folder (e.g.: `E:\210527_NDX550487_RUO_0008_AHTWG2AFX2\Alignment_1\20210528_061936`).
  It contains the two files `CompletedJobInfo.xml` and `SampleSheetUsed.csv` and a subdirectory named `Fastq` containing all the`.fastq.gz` sequencing files.
- `--outdir` : directory where to create the output

This command is the **first step** (analyzing `CompletedJobInfo.xml` and `SampleSheetUsed.csv` ):

- It will generate a samples TSV file in the output directory: `working/samples/samples.`_{some date}_`.tsv`.
  You can copy the content of this file into your `samples.tsv`
- It will generate a file `working/samples/movedata.sh`, you can run it with:
  ```console
  bash working/samples/movedata.sh
  ```
  and that file will in turn perform the **second step**: hard-linking all the files from `downloads/20210528_061936/Fastq` into `working/samples/`_{sample name}_`/`_{some date}_.
- once the second step has been performed, you should be able to run V-pipe.
- if the option `--batch` is provided, it will also generate an additional file `working/samples/batch.`_{some date}_`.yaml` including extra information gathered from the files (e.g. the _Library preparation kit_ listed in the input CSV). The parameter of `--batch` is used to provide the name of the lab for the `lab:` field in this file.

Available options:

```console
./sort_samples_jobinfo  -h
usage: sort_samples_jobinfo [-h] -S DIR [-f DIR] [-o DIR] [-m MODE] [-L CPLINK] [-b LAB] [-s] [-a] [-l] [-p TSV]

Uses CompletedJobInfo.xml and SampleSheetUsed.csv from Illumina Analysis Software

optional arguments:
  -h, --help            show this help message and exit
  -S DIR, --sourcedir DIR
                        directory containing CompletedJobInfo.xml and SampleSheetUsed.csv
  -f DIR, --fastqdir DIR
                        directory containing .fastq.gz files if not in 'Fastq' subdirectory
  -o DIR, --outdir DIR  output directory
  -m MODE, --mode MODE  POSIX file access mode to be passed to mkdir
  -L CPLINK, --linking CPLINK
                        parameter to pass to `cp` for linking files instead of copying their data
  --force               Force overwriting any existing file when moving
  -b LAB, --batch LAB   generate batch description
  -s, --summary         Only display a summary of datasets, not an exhaustive list of all samples
  -a, --append          Append to the end of movedatafiles.sh, instead of overwritting (use when calling from an external combiner wrapper)
  -g, --staging         Write samples list in .tsv.staging and only rename them to the final .tsv at the end of movedatafiles.sh if there were no errors.
  -l, --forcelanes      Explicitly look for sample in each lane (for replicates across lanes)
  -p TSV, --patchmap TSV
                        patchmap file to rename samples
```
