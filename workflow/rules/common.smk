__author__ = "Susana Posada-Cespedes"
__author__ = "David Seifert"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"

import csv
import os
import pathlib
import re
import typing

from datetime import datetime, timezone

# for legacy
import configparser
import json  # NOTE jsonschema  is only part of the full snakemake, not snakemake-minimal

from collections import UserDict

from snakemake.io import load_configfile
from snakemake.utils import update_config, validate, min_version

from snakemake.sourcecache import infer_source_file
from snakemake.common import is_local_file

import logging

LOGGER = logging.getLogger("snakemake.logging")

if not "VPIPE_BENCH" in dir():
    VPIPE_BENCH = False

##########################
#   Configuration file   #
##########################

# even if the --configfile option overload is used this default *must* always exist:
# - the _content_ of these extra configfiles will overwrite the basefile's _content_
# - not the _filenames_ themselves
if os.path.exists("config/config.yaml"):

    configfile: "config/config.yaml"

elif os.path.exists("config.yaml"):

    configfile: "config.yaml"


# NOTE the remote URL relies on snakemake sourcecache features (e.g.: infer_source_file) introduced in 6.8.1
# NOTE runtime resources calculation introduced in 7.11.0
min_version("7.11.0")


def cacheopen(fname, mode="r", localsource=False):
    """
    get a file handle of the fname.
    snakemake.sourcecache will automatically cache it on-the-fly
    """

    return open(cachepath(fname, executable=False, localsource=localsource), mode)


def cachepath(fname, executable=False, localsource=False):
    """
    get the local filename for a given file
     - local files: well it's just the filename it-self
     - remote URLs: they get cached by snakemake.sourcecache
           and we use the cache's filename
    """

    # TODO revisit official snakemake documented way to access remote source files when it maturees.

    # yes, srcdir duplicates work done later by workflow.source_path, but it makes is_local_file run on the full prefixed path.
    fullpath = srcdir(fname) if localsource else fname
    if is_local_file(fullpath):
        # local files as-is
        return fullpath

    # remote: cache (if not present yet)
    cached = (
        workflow.source_path(fname)
        if localsource
        else workflow.sourcecache.get_path(infer_source_file(fullpath))
    )

    if executable:
        # r to x (e.g.: u+rw,g+r,o-rwx 0o640 => 0o750 u+rwx,g+rx,o-rwx)
        mode = os.stat(cached).st_mode
        mode |= (mode & 0o444) >> 2
        os.chmod(cached, mode)

    # HACK currently, sourcecache doesn't handle dates and cached file often are from "today". This trick avoid constantly re-running rules with http-fetched dependencies
    return ancient(cached)


def load_legacy_ini(ininame, schemaname):
    """
    Load a legacy INI-style config file with Python configparser
    and uses a schema to load it correctly:

    - configparser does not guess the datatpye of value in config files
      always storing them internally as strings.
      it provides getter methods for integers, floats and booleans.

      we apply the getter corresponding to the type expected in the schema

    - configpasrer keys are not case-sensitive and stored in lowercase,
      whereas JSON and YAML are case sensitive.

      we convert to the case expected by the schema
    """
    cp = configparser.ConfigParser()
    cp.read(ininame)

    with cacheopen(schemaname, mode="rt") as sf:
        sch = json.load(sf)

    ini = {}
    # build the result structure by:
    # - looping the whole configparser object
    # - applying types from the jsonschema
    # - fixing the case of keys
    for name, section in cp.items():
        if name == "DEFAULT":
            # configparser-specific, skip
            continue
        if name not in sch["properties"]:
            # technically, missing part should be handled by validator,
            # but we need valid names anyway to do the types
            raise ValueError(
                "Invalid section in {inif} :".format(sec=name, inif=ininame), name
            )
        ini[name] = {}
        # configpasrer keys are not case-sensitive and stored in lowercase
        lcentries = {k.lower(): k for k in sch["properties"][name]["properties"].keys()}
        for entry, value in section.items():
            if entry not in lcentries:
                raise ValueError(
                    "Invalid entry in section {sec} of {inif} :".format(
                        sec=name, inif=ininame
                    ),
                    entry,
                )
            entry = lcentries[entry]
            # configparser doesn't convert type automatically
            # convert if necessary based on schema
            ty = sch["properties"][name]["properties"][entry]["type"]
            if ty == "boolean":
                value = section.getboolean(entry)
            elif ty == "integer":
                value = section.getint(entry)
            elif ty == "number":
                value = section.getfloat(entry)
            # string are left as-is
            ini[name][entry] = value
    return ini


def process_config(config):
    """
    process configuration.

    - precedence logic:
      snakemake (configfile(s) + --config) >> legacy configparse INI >> virus base config >> schema default

    - validate with schema
    """

    schema = srcdir("../schemas/config_schema.json")

    # merging of legacy INI-style
    if os.path.exists("vpipe.config"):
        LOGGER.info("Importing legacy configuration file vpipe.config")
        # snakemake configuration overwrites legacy vpipe.config
        cur_config = config
        config = load_legacy_ini("vpipe.config", schema)
        update_config(config, cur_config)

    # merging of virus' base configuration
    vf = None
    try:
        # search for file location
        if config["general"]["virus_base_config"]:
            # shorthand - e.g.: hiv
            try:
                vf = cacheopen(
                    "../../config/{VIRUS}.yaml".format(
                        VIRUS=config["general"]["virus_base_config"],
                    ),
                    mode="rt",
                    localsource=True,
                )
            except FileNotFoundError:
                vf = None

            # normal search (with normal macro expansion, like the default values)
            if vf is None:
                try:
                    vf = cacheopen(
                        config["general"]["virus_base_config"].format(
                            VPIPE_BASEDIR=VPIPE_BASEDIR
                        ),
                        mode="rt",
                    )
                except FileNotFoundError:
                    vf = None

            if vf is None:
                raise ValueError(
                    "Cannot find virus base config",
                    config["general"]["virus_base_config"],
                )
    # (empty entry)
    except TypeError:
        vf = None
    except KeyError:
        vf = None

    if vf is not None:
        # current configuration overwrites virus base config
        cur_config = config
        config = load_configfile(vf)
        vf.close()
        if "name" in config:
            LOGGER.info("Using base configuration virus %s" % config.pop("name"))
        else:
            LOGGER.info(
                "Using base configuration from %s"
                % cur_config["general"]["virus_base_config"]
            )
        update_config(config, cur_config)
    else:
        LOGGER.info("No virus base configuration, using defaults")

    # validates, but also fills up default values:
    validate(config, schema, set_default=True)
    # use general.threads entry as default for all affected sections
    # if not specified:
    for name, section in config.items():
        if name == "general":
            continue
        if not isinstance(section, dict):
            continue
        if "threads" not in section:
            section["threads"] = config["general"]["threads"]

        # interpolate some parameters
        # (currently only the base directory for resources packaged in V-pipe)
        for entry, value in section.items():
            if isinstance(value, str):
                section[entry] = value.format(VPIPE_BASEDIR=VPIPE_BASEDIR)

    """
    The following hack supports `.` access to recursive dictionary entries,
    e.g.  config.input instead of config["input"] so we don't have to change
    all existing rules.
    """

    def wrap(dd):
        if isinstance(dd, dict):
            dd = UserDict(dd)
            for k, v in dd.items():
                dd[k] = wrap(v)
            dd.__dict__.update(dd)
        return dd

    return wrap(config)


config = process_config(config)


def load_protocols(pyaml):
    """
    Load the look-up YAML file specifying per-protocol specific, usefull for a 4th column in the samples TSV table
    """
    if not pyaml:
        return {}

    # normal search (with normal macro expansion, like the default values)
    try:
        pf = cacheopen(
            pyaml.format(VPIPE_BASEDIR=VPIPE_BASEDIR),
            mode="rt",
        )
    except FileNotFoundError as e:
        raise ValueError(
            f"Cannot find protocols look-up YAML file {pyaml}",
        ) from e

    py = load_configfile(pf)

    # interpolate some parameters
    # (currently only the base directory for resources packaged in V-pipe)
    for name, section in py.items():
        if not isinstance(section, dict):
            continue
        for entry, value in section.items():
            if isinstance(value, str):
                section[entry] = value.format(VPIPE_BASEDIR=VPIPE_BASEDIR)

    return py


protocols = load_protocols(config["input"]["protocols_file"])


def cohortdir(fname):
    return os.path.normpath(
        os.path.join(
            config.output["datadir"],
            config.output["cohortdir"],
            fname,
        )
    )


# handle Genes GFF
if config.input["gff_directory"]:
    assert is_local_file(
        config.input["gff_directory"]
    ), f"""ERROR: section 'input' property 'gff_directory' cannot be a remote URL {config.input['gff_directory']}.
when not running V-pipe locally (e.g. snakedeploy) you need to either:
- set gff_directory to empty ''
- copy to a local directory and specify that as the gff_directory"""

if (
    "genes_gff" in config.frameshift_deletions_checks
    and config.frameshift_deletions_checks["genes_gff"]
):
    # HACK legacy behaviour
    config.input["genes_gff"] = config.frameshift_deletions_checks["genes_gff"]
    LOGGER.warning(
        f"WARNING: Genes .gff file {config.input['genes_gff']} specified in section 'frameshift_deletions_checks'. This is a legacy option and will be deprecated in the future. Please move property to section 'input'."
    )
elif config.input["genes_gff"]:
    LOGGER.info("Genes are specified using file %s" % config.input["genes_gff"])

# HACK not specified, trying to pick one.
if not config.input["genes_gff"] and config.input["gff_directory"]:
    gtf_glob = ".[Gg][FfTt][Ff]*"
    LOGGER.warning(
        "WARNING: No genes .gff file specified in section 'input' property 'genes_gff'."
    )
    import glob

    mf = (
        load_protocols(config.input["metainfo_file"])
        if config.input["metainfo_file"]
        else None
    )

    if mf and "gff" in mf:
        # HACK creative abuse of protocol loader
        LOGGER.warning(f"Searching {config.input['metainfo_file']} metainfo_file...")

        # Search for "gene" in the desciption
        for fn, dsc in mf["gff"].items():
            if "gene" in dsc.lower():
                srch = glob.glob(
                    os.path.join(config.input["gff_directory"], f"{fn}{gtf_glob}")
                )
                if len(srch) >= 1:
                    config.input["genes_gff"] = srch[0]
                    LOGGER.warning(
                        f"WARNING: Autoselecting {config.input['genes_gff']} based on desciption <{dsc}>"
                    )
                    break

    if not config.input["genes_gff"]:
        LOGGER.warning(f"Searching gff_directory {config.input['gff_directory']}...")
        srch = glob.glob(
            os.path.join(config.input["gff_directory"], f"*[Gg]ene*{gtf_glob}")
        )
        if len(srch) >= 1:
            srch.sort()
            config.input["genes_gff"] = srch[0]
            LOGGER.warning(
                f"WARNING: Autoselecting {config.input['genes_gff']} based on file name."
            )
        elif mf and "gff" in mf:
            fn, dsc = list(mf["gff"].items())[0]
            srch = glob.glob(os.path.join(config.input["gff_directory"], f"{fn}{gtf_glob}"))
            if len(srch) >= 1:
                srch.sort()
                config.input["genes_gff"] = srch[0]
                LOGGER.warning(
                    f"WARNING: Autoselecting {config.input['genes_gff']} <{dsc}> from directory."
                )

    if not config.input["genes_gff"]:
        srch = glob.glob(os.path.join(config.input["gff_directory"], f"*{gtf_glob}"))
        if len(srch) >= 1:
            srch.sort()
            config.input["genes_gff"] = srch[0]
            LOGGER.warning(
                f"WARNING: Autoselecting {config.input['genes_gff']} from directory."
            )

if not config.input["genes_gff"]:
    LOGGER.warning(
        f"WARNING: No genes .gff file specified in section 'input' property 'genes_gff'. Some reports will be degraded. Please provide one if possible to improve reporting."
    )


########################
#   Samples TSV file   #
########################

# 2. glob patients/samples + store as TSV if file is not provided

# if file containing samples exists, proceed
# to build list of target files
if not os.path.isfile(config.input["samples_file"]):
    # sample file does not exist, have to first glob
    # all samples' data and then construct sample
    # list that would pass QA checks
    LOGGER.info(
        f"Missing Sample list file {config.input['samples_file']}, automatically scanning {config.input['datadir']} for samples..."
    )

    sample_pairs = glob_wildcards(
        os.path.join(config.input["datadir"], "{sample_date}", "raw_data", "{file}")
    )

    if not set(sample_pairs.sample_date):
        LOGGER.warning(
            f"WARNING: No samples found in {config.input['datadir']}. Not generating {config.input['samples_file']}."
        )
    else:
        c = 0
        with open(config.input["samples_file"], "w") as outfile:
            for i in set(sample_pairs.sample_date):
                (sample, date) = (
                    [x.strip() for x in i.split(os.path.sep) if x.strip()] + [""]
                )[:2]
                outfile.write("{}\t{}\n".format(sample, date))
                c += 1
        LOGGER.info(f"{c} samples found in {config.input['datadir']}.")


# TODO: have to preprocess samples files to filter likely failures
# 1.) Determine 5%/95% length of FASTQ files
# 2.) Determine whether FASTQ would survive
# 3.) Detect amplicon protocols using ampseer
# NOTE: Currently, "1.)" is partially handled in utils/sort_sample_dumb.

# 3. load samples from TSV and create list of samples
#
# This list is reused on subsequent runs

sample_list = []
sample_table = {}
sample_record = typing.NamedTuple("sample_record", [("sample_id", str), ("date", str)])
sample_id_patchmap = {}
sample_dir = {}  # directory => samples record
sample_paths = {}  # sample record => dir
sample_proto_count = 0
sample_default_count = 0
sample_1level_count = 0
sample_2level_count = 0
sample_1level_line = 0
sample_2level_line = 0


def guess_sample(path):
    """
    function to guess the current sample's record from the path
    e.g.:
      "results/Sam12/seq_20211108/"
       =>
      samplet_id="Sam12", date="seq_20211108"
    """
    s_rec = sample_record(
        **dict(
            zip(
                ["sample_id", "date"],
                list(os.path.normpath(path).split(os.path.sep)[-2:]),
            )
        )
    )
    # HACK to handle gracefully non two-level samples (e.g.: single level)
    # TODO replace with string suffix search from sample_dir
    return sample_id_patchmap.get(s_rec, s_rec)


# TODO this should in the end go into some structure, ideally a peppy (Python module for Portable Ecapsulated Project)
sample_row = typing.NamedTuple(
    "sample_row", [("num", int), ("len", int), ("protocol", str)]
)


if not os.path.isfile(config.input["samples_file"]):
    LOGGER.warning(f"WARNING: Sample list file {config.input['samples_file']} not found.")
else:
    with open(config.input["samples_file"], newline="") as csvfile:
        spamreader = csv.reader(csvfile, delimiter="\t")

        for row in spamreader:
            if len(row) == 0 or row[0][0] == "#":
                # Skip completely empty lines, or comments begining with '#' (numpy style)
                continue
            assert (
                len(row) >= 2
            ), "ERROR: Line '{}' does not contain at least two entries!".format(
                spamreader.line_num
            )
            sample_tuple = sample_record(sample_id=row[0], date=row[1])
            sample_list.append(sample_tuple)

            assert (
                sample_tuple not in sample_table
            ), "ERROR: on line {} sample '{}-{}' is not unique, previous was on line {}".format(
                spamreader.line_num, row[0], row[1], sample_table[sample_tuple].num
            )

            if not (sample_tuple.sample_id and sample_tuple.date):
                if not sample_1level_count:
                    sample_1level_line = spamreader.line_num
                sample_1level_count += 1
                # HACK to handle gracefully non two-level samples (e.g.: single level)

                # guess_sample can return wrong guesses if there is only 1 level: the output directory will be assigned to sample_id
                patch_tuple = guess_sample(
                    os.path.join(
                        config.output["datadir"],
                        sample_tuple.sample_id,
                        sample_tuple.date,
                    )
                )
                # we keep track of such mis-haps so we can patch them.
                sample_id_patchmap[patch_tuple] = sample_tuple
            else:
                if not sample_2level_count:
                    sample_2level_line = spamreader.line_num
                sample_2level_count += 1

            # defaults (if columns are missing in TSV)
            l = config.input[
                "read_length"
            ]  # All samples are assumed to have same read length and the default, 250 bp
            p = None  # protocol-specific are assumed to be passed by config options

            if (len(row) >= 3) and row[2]:
                # Extract read length from input.samples_file. Samples may have
                # different read lengths. Reads will be filtered out if read length
                # after trimming is less than trim_cutoff * read_length.
                try:
                    l = int(row[2])
                except ValueError as e:
                    raise ValueError(
                        "ERROR: Wrong read-length value given on line {} for sample '{}-{}'. If present, third column in a TSV file MUST be a number. Not <{}>. Please read: config/README.md or https://github.com/cbg-ethz/V-pipe/tree/master/config".format(
                            spamreader.line_num, row[0], row[1], row[2]
                        )
                    ) from e

            if (len(row) >= 4) and row[3]:
                # Extract protocol name from sample file.
                # Over the time of a long-running experiment, protocols can change to adapt to changing mix of present variants
                # e.g.: Primers might change due to SNVs in new variants
                p = row[3]

                if not len(protocols):
                    raise ValueError(
                        "ERROR: Protocol short name <{}> specified on line {} for sample '{}-{}', but no protocols_file defined in section 'input' of the configuration. If present, fourth column in a TSV file MUST be a shortname specified in the protocols YAML look-up file. Please read: config/README.md or https://github.com/cbg-ethz/V-pipe/tree/master/config".format(
                            p, spamreader.line_num, row[0], row[1]
                        )
                    )

                if p not in protocols:
                    raise ValueError(
                        "ERROR: Wrong protocol short name <{}> on line {} for sample '{}-{}'. If present, fourth column in a TSV file MUST be a shortname specified in the protocols look-up file: [{}]. Please read: config/README.md or https://github.com/cbg-ethz/V-pipe/tree/master/config".format(
                            p,
                            spamreader.line_num,
                            row[0],
                            row[1],
                            (";".join(protocols.keys())),
                        )
                    )

                sample_proto_count += 1
            else:
                sample_default_count += 1

            sample_table[sample_tuple] = sample_row(
                num=spamreader.line_num, len=l, protocol=p
            )

if config["output"]["trim_primers"]:
    if sample_default_count and not config["input"]["primers_bedfile"]:
        raise ValueError(
            "ERROR: {} sample(s) do not specify any protocol short name. If there is no fourth column in a TSV file, a default primers_bedfile MUST be specified in the config file in section input, option primers_bedfile. Please read: config/README.md or https://github.com/cbg-ethz/V-pipe/tree/master/config".format(
                sample_default_count
            )
        )
    elif (0 == sample_default_count) and config["input"]["primers_bedfile"]:
        LOGGER.warning(
            "NOTE: no sample uses the default primers_bedfile, each sample specifies a protocol short name."
        )

if len(protocols) and (0 == sample_proto_count):
    LOGGER.warning(
        "WARNING: protocols YAML look-up file <{}> specified, but no sample ever uses it: fourth column absent from samples TSV file.".format(
            config["input"]["protocols_file"]
        )
    )

if sample_1level_count and sample_2level_count:
    LOGGER.warning(
        "WARNING: samples TSV contains both {} samples with 2-level hierarchy (starting at line {}) and {} samples with 1-level hierarchy -- i.e. one of the two columns is empty (starting at line {}). Such mixing isn't thoroughly tested, consider it unsupported. All bugs encountered are features ;-)".format(
            sample_2level_count, sample_2level_line, sample_1level_count, sample_1level_line
        )
    )


# 4. generate list of target files
all_files = []
alignments = []
vicuna_refs = []
references = []
consensus = []
trimmed_files = []
fastqc_files = []
results = []
visualizations = []
diversity_measures = []
datasets = []
IDs = []
dehumanized_raw_reads = []
upload_markers = []

# do we use raw alignments or trimmed alignments?
alignment_file = (
    "alignments/REF_aln.bam"
    if not config.output["trim_primers"]
    else "alignments/REF_aln_trim.bam"
)
alignment_wildcard = "{dataset}/" + alignment_file

for srec in sample_list:
    # WARNING the following makes sure to gracefully handle trailing slashes in the user-provided paths in datadir
    sdir = os.path.join(
        config.output["datadir"], os.path.normpath(os.path.join(srec.sample_id, srec.date))
    )
    sample_dir[sdir] = srec
    sample_paths[srec] = sdir

    alignments.append(os.path.join(sdir, alignment_file))
    # if config.output["QA"]:
    #    alignments.append(os.path.join(sdir, "QA_alignments/coverage_ambig.tsv"))
    #    alignments.append(os.path.join(sdir, "QA_alignments/coverage_majority.tsv"))

    vicuna_refs.append(os.path.join(sdir, "references/vicuna_consensus.fasta"))
    references.append(os.path.join(sdir, "references/ref_"))

    consensus.append(os.path.join(sdir, "references/ref_ambig.fasta"))
    consensus.append(os.path.join(sdir, "references/ref_ambig_dels.fasta"))
    consensus.append(os.path.join(sdir, "references/ref_majority.fasta"))
    consensus.append(os.path.join(sdir, "references/ref_majority_dels.fasta"))

    consensus.append(os.path.join(sdir, "references/consensus.bcftools.fasta"))

    if config.output["QA"]:
        alignments.append(os.path.join(sdir, "references/ref_majority_dels.matcher"))
        alignments.append(os.path.join(sdir, "references/frameshift_deletions_check.tsv"))

    trimmed_files.append(os.path.join(sdir, "preprocessed_data/R1.fastq.gz"))
    if config.input["paired"]:
        trimmed_files.append(os.path.join(sdir, "preprocessed_data/R2.fastq.gz"))

    fastqc_files.append(os.path.join(sdir, "extracted_data/R1_fastqc.html"))
    if config.input["paired"]:
        fastqc_files.append(os.path.join(sdir, "extracted_data/R2_fastqc.html"))

    datasets.append(sdir)
    IDs.append(
        config.general["id_separator"].join([srec.sample_id, srec.date])
        if srec.sample_id and srec.date
        else srec.sample_id or srec.date
    )

    # SNV
    if config.output["snv"]:
        # in adition to standard VCF files, ShoRAH2 also produces CSV tables
        if config.general["snv_caller"] == "shorah":
            results.append(os.path.join(sdir, "variants/SNVs/snvs.csv"))
        elif config.general["snv_caller"] == "viloca":
            results.append(
                os.path.join(sdir, "variants/SNVs/snv/cooccurring_mutations.csv")
            )
        # all snv callers ('shorah', 'lofreq') produce standard VCF files
        results.append(os.path.join(sdir, "variants/SNVs/snvs.vcf"))
    # local haplotypes
    if config.output["local"]:
        if config.general["snv_caller"] == "shorah":
            results.append(os.path.join(sdir, "variants/SNVs/snvs.csv"))
        elif config.general["snv_caller"] == "viloca":
            results.append(
                os.path.join(sdir, "variants/SNVs/snv/cooccurring_mutations.csv")
            )
    # global haplotypes
    if config.output["global"]:
        if config.general["haplotype_reconstruction"] == "savage":
            results.append(os.path.join(sdir, "variants/global/contigs_stage_c.fasta"))
        elif config.general["haplotype_reconstruction"] == "haploclique":
            results.append(os.path.join(sdir, "variants/global/quasispecies.bam"))
        elif config.general["haplotype_reconstruction"] == "predicthaplo":
            if config.input["paired"]:
                results.append(
                    os.path.join(sdir, "variants/global/predicthaplo_haplotypes.fasta")
                )
            else:
                raise NotImplementedError("PredictHaplo only works with paired-end reads")

    # visualization
    if not config.output["snv"] and config.output["visualization"]:
        raise RuntimeError(
            "Cannot generate visualization without calling variants (make sure to set `snv = True`) in config."
        )

    if config.output["snv"] and config.output["visualization"]:
        visualizations.append(os.path.join(sdir, "visualization/snv_calling.html"))
        visualizations.append(os.path.join(sdir, "visualization/alignment.html"))

    # upload related stuff
    if config.output["dehumanized_raw_reads"]:
        dehumanized_raw_reads.append(os.path.join(sdir, "raw_uploads", "dehuman.cram"))

    if config.output["upload"]:
        upload_markers.append(os.path.join(sdir, "upload_prepared.touch"))

# merge lists containing expected output
all_files = (
    alignments
    + consensus
    + results
    + visualizations
    + dehumanized_raw_reads
    + upload_markers
)

# diversity measures
if not config.output["snv"] and config.output["diversity"]:
    raise RuntimeError(
        "Cannot generate diversity without calling variants (make sure to set `snv = True`) in config."
    )

if config.output["snv"] and config.output["diversity"]:
    all_files.append(cohortdir("aggregated_diversity.csv"))


IDs = ",".join(IDs)


# 5. Locate reference and parse reference identifier
def get_reference_name(reference_file):
    with open(reference_file, "r") as infile:
        reference_name = infile.readline().rstrip()
    reference_name = reference_name.split(">")[1]
    reference_name = reference_name.split(" ")[0]
    return reference_name


if not VPIPE_BENCH:
    reference_file = config["input"]["reference"]
    if not reference_file:
        raise ValueError(
            f"ERROR: No input reference in configuration. Please read: config/README.md or https://github.com/cbg-ethz/V-pipe/tree/master/config"
        )
    elif not is_local_file(reference_file):
        reference_file_alt = cachepath(reference_file)
        LOGGER.info(f"Caching {reference_file} into {reference_file_alt}")
        reference_file = reference_file_alt
        reference_name = get_reference_name(reference_file_alt)
    elif not os.path.isfile(reference_file):
        reference_file_alt = os.path.join("references", reference_file)
        LOGGER.warning(
            f"WARNING: Reference file {reference_file} not found. Trying {reference_file_alt}."
        )
        reference_file = reference_file_alt
        if not os.path.isfile(reference_file):
            raise ValueError(f"ERROR: Reference file {reference_file} not found.")
        reference_name = get_reference_name(reference_file)
    else:
        reference_name = get_reference_name(reference_file)

    if config.output["dehumanized_raw_reads"] and not os.path.isfile(
        config.dehuman["ref_host"]
    ):
        LOGGER.warning(
            f"WARNING: Host organism reference file {config.dehuman['ref_host']} not found, it will be downloaded from {config.dehuman['ref_host_url']}."
        )


# Auxiliary functions

# TODO These shoudle eventually go into additional columns once we move to proper dataset


def protocol_proto_option(proto, option, s_rec=None):
    """return specific option given the provided proto"""
    # no sample-specific protocol => use from config
    if not proto:
        return config.input[option]

    #  samples-sepcific protocol => pick from the protocol look-up YAML
    try:
        return protocols[proto][option]
    except (KeyError, TypeError) as e:
        raise KeyError(
            f"no {option} defined for protocol {proto} used by sample {s_rec.sample_id}-{s_rec.date}"
            if s_rec
            else f"no {option} defined for protocol {proto}"
        ) from e


def get_sample_protocol(wildcards):
    """get the 4th column (protocol) of the current sample"""
    s_rec = guess_sample(wildcards.dataset)
    proto = sample_table[s_rec].protocol

    return proto


def protocol_option(wildcards, option):
    """get a per-protocol specific configuration option, depending on the current sample's 4th column"""
    # skip if no sample ever has 4th column
    if 0 == sample_proto_count:
        return config.input[option]

    proto = get_sample_protocol(wildcards)

    return protocol_proto_option(proto, option)


def ID(wildcards):
    s_rec = guess_sample(wildcards.dataset)
    try:
        # normal two-level
        return (
            config.general["id_separator"].join(s_rec)
            if s_rec.sample_id and s_rec.date
            else s_rec.sample_id or s_rec.date
        )
    except TypeError:
        # HACK single-level
        return s_rec.sample_id or s_rec.date


def window_lengths(wildcards):
    window_len = []
    for s in sample_list:
        read_len = sample_table[s].len
        aux = int((read_len * 4 / 5 + config.snv["shift"]) / config.snv["shift"])
        window_len.append(str(aux * config.snv["shift"]))

    return ",".join(window_len)


def shifts(wildcards):
    shifts = []
    for s in sample_list:
        read_len = sample_table[s].len
        aux = int((read_len * 4 / 5 + config.snv["shift"]) / config.snv["shift"])
        shifts.append(str(aux))

    return ",".join(shifts)


def get_maxins(wildcards):
    if "maxins" in config["bowtie_align"]:
        return config.bowtie_align["maxins"]
    else:
        s_rec = guess_sample(wildcards.dataset)
        read_len = sample_table[s_rec].len
        return 4 * read_len


# WARNING needs to handle trailing slashes gracefully and re-use the exact two-levels of the directory structure (patient / date)
def rebase_datadir(base, dataset):
    s_rec = guess_sample(dataset)
    try:
        # normal two-level
        return os.path.join(base, *list(s_rec))
    except TypeError:
        # HACK single-level
        return os.path.join(base, s_rec.sample_id or s_rec.date)


def raw_data_file(wildcards, pair):
    indir = os.path.join(
        rebase_datadir(config.input["datadir"], wildcards.dataset), "raw_data"
    )
    list_output = []
    rx_fq = re.compile(
        r"[^/]*{}\.(fastq\.gz|fastq|fq|fq\.gz)$".format(
            ("R%u%s" % (pair, config.input["fastq_suffix"])) if pair else ""
        )
    )
    # if directory was deleted: no input files to restart alignement
    if not os.path.exists(indir):
        return []
    # otherwise, scan for .fastq files
    for p in os.listdir(indir):
        if rx_fq.search(p):
            list_output.append(os.path.join(indir, p))

    if len(list_output) == 0:
        raise ValueError(
            "Missing input files for sample in: {} - Unexpected file name?".format(indir)
        )

    return list_output


# TODO replace with raw_data_file
def construct_input_fastq(wildcards):
    indir = os.path.join(
        rebase_datadir(config.input["datadir"], wildcards.dataset), "raw_data"
    )
    aux = glob_wildcards(
        indir + "/{prefix, [^/]+}" + "{ext, (\.fastq|\.fastq\.gz|\.fq|\.fq\.gz)}"
    )
    if config.input["paired"]:
        inferred_values = glob_wildcards(
            indir + "/{file}R" + wildcards.pair + config.input["fastq_suffix"] + aux.ext[0]
        )
    else:
        inferred_values = glob_wildcards(indir + "/{file}" + aux.ext[0])

    list_output = []
    file_extension = aux.ext[0].split(".gz")[0]
    for i in inferred_values.file:
        if config.input["paired"]:
            list_output.append(
                os.path.join(
                    indir,
                    "".join(
                        (
                            i,
                            "R",
                            wildcards.pair,
                            config.input["fastq_suffix"],
                            file_extension,
                        )
                    ),
                )
            )
        else:
            list_output.append(os.path.join(indir, "".join((i, file_extension))))
    if len(list_output) == 0:
        raise ValueError(
            "Missing input files for rule extract: {}/raw_data/ - Unexpected file name?".format(
                wildcards.dataset
            )
        )

    return list_output


def temp_prefix(p):
    return os.path.join(config.general["temp_prefix"], p)


def temp_with_prefix(p):
    return temp(temp_prefix(p))


###############
#   Globals   #
###############

functions = cachepath(
    "../scripts/functions.sh",
    executable=True,
    localsource=True,
)
