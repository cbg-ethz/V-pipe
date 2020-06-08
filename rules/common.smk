import csv
import collections
import configparser
import os
import typing

__author__ = "Susana Posada-Cespedes"
__author__ = "David Seifert"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"

VPIPE_BENCH = True if os.environ.get('VPIPE_BENCH') is not None else False


# 1. parse config file
#
# the config object sets defaults for a number of parameters and
# validates parameters given by the user to be sensible

VPIPE_CONFIG_OPTS = True if os.environ.get('VPIPE_CONFIG_OPTS') is not None else False

class VpipeConfig(object):
    'Class used to encapsulate the configuration properties used by V-pipe'

    __RECORD__ = typing.NamedTuple(
        "__RECORD__", [('value', typing.Any), ('type', type)])
    __MEMBER_DEFAULT__ = collections.OrderedDict([
        ('general', {
            'threads': __RECORD__(value=4, type=int),
            'aligner': __RECORD__(value='ngshmmalign', type=str),
            'snv_caller': __RECORD__(value='shorah', type=str),
            'haplotype_reconstruction': __RECORD__(value='savage', type=str)
        }),
        ('input', {
            'datadir': __RECORD__(value='samples', type=str),
            'samples_file': __RECORD__(value='samples.tsv', type=str),
            'paired': __RECORD__(value=True, type=bool),
            'fastq_suffix': __RECORD__(value='', type=str),
            'trim_percent_cutoff': __RECORD__(value=0.8, type=float),
            'reference': __RECORD__(value='references/HXB2.fasta', type=str),
            'gff_directory': __RECORD__(value='', type=str),
            'primers_file': __RECORD__(value='', type=str),
        }),
        ('output', {
            'QA': __RECORD__(value=False, type=bool),
            'snv': __RECORD__(value=True, type=bool),
            'local': __RECORD__(value=True, type=bool),
            'global': __RECORD__(value=True, type=bool),
            'visualization': __RECORD__(value=True, type=bool),
        }),
        ('applications', {
            'gunzip': __RECORD__(value="gunzip", type=str),
            'prinseq': __RECORD__(value="prinseq-lite.pl", type=str),
            'fastqc': __RECORD__(value="fastqc", type=str),
            'vicuna': __RECORD__(value="vicuna", type=str),
            'indelfixer': __RECORD__(value="InDelFixer", type=str),
            'consensusfixer': __RECORD__(value="ConsensusFixer", type=str),
            'picard': __RECORD__(value="picard", type=str),
            'bwa': __RECORD__(value="bwa", type=str),
            'bowtie_idx': __RECORD__(value="bowtie2-build", type=str),
            'bowtie': __RECORD__(value="bowtie2", type=str),
            'samtools': __RECORD__(value="samtools", type=str),
            'extract_consensus': __RECORD__(value="extract_consensus", type=str),
            'mafft': __RECORD__(value="mafft", type=str),
            'ngshmmalign': __RECORD__(value="ngshmmalign", type=str),
            'convert_reference': __RECORD__(value="convert_reference", type=str),
            'extract_seq': __RECORD__(value="extract_seq", type=str),
            'coverage_stats': __RECORD__(value="coverage_stats", type=str),
            'remove_gaps_msa': __RECORD__(value="remove_gaps_msa", type=str),
            'minority_freq': __RECORD__(value="minority_freq", type=str),
            'extract_coverage_intervals': __RECORD__(value="extract_coverage_intervals", type=str),
            'shorah': __RECORD__(value="shorah shotgun", type='str'),
            'lofreq': __RECORD__(value="lofreq", type=str),
            'bcftools': __RECORD__(value="bcftools", type=str),
            'haploclique': __RECORD__(value="haploclique", type='str'),
            'compute_mds': __RECORD__(value="compute_mds", type='str'),
            'savage': __RECORD__(value="savage", type='str')
        }),

        ('gunzip', {
            'mem': __RECORD__(value=30000, type=int),
            'time': __RECORD__(value=60, type=int),
        }),
        ('extract', {
            'mem': __RECORD__(value=10000, type=int),
            'time': __RECORD__(value=20, type=int),
        }),
        ('preprocessing', {
            'mem': __RECORD__(value=2000, type=int),
            'time': __RECORD__(value=235, type=int),
            'conda': __RECORD__(value='', type=str),

            'extra': __RECORD__(value='-ns_max_n 4 -min_qual_mean 30 -trim_qual_left 30 -trim_qual_right 30 -trim_qual_window 10', type=str),
        }),
        ('fastqc', {
            'mem': __RECORD__(value=2000, type=int),
            'time': __RECORD__(value=235, type=int),
            'conda': __RECORD__(value='', type=str),
            'threads': __RECORD__(value=0, type=int),

            'no_group': __RECORD__(value=False, type=bool),
        }),
        ('initial_vicuna', {
            'mem': __RECORD__(value=1000, type=int),
            'time': __RECORD__(value=600, type=int),
            'threads': __RECORD__(value=0, type=int),
            'conda': __RECORD__(value='', type=str),
        }),
        ('initial_vicuna_msa', {
            'mem': __RECORD__(value=10000, type=int),
            'time': __RECORD__(value=235, type=int),
            'threads': __RECORD__(value=0, type=int),
            'conda': __RECORD__(value='', type=str),
        }),
        ('create_vicuna_initial', {
            'conda': __RECORD__(value=f'{VPIPE_BASEDIR}/envs/smallgenomeutilities.yaml', type=str),
        }),
        ('hmm_align', {
            'mem': __RECORD__(value=1250, type=int),
            'time': __RECORD__(value=1435, type=int),
            'threads': __RECORD__(value=0, type=int),
            'conda': __RECORD__(value='', type=str),

            'leave_msa_temp': __RECORD__(value=False, type=bool),
            'extra': __RECORD__(value='', type=str),
        }),
        ('sam2bam', {
            'mem': __RECORD__(value=5000, type=int),
            'time': __RECORD__(value=30, type=int),
            'conda': __RECORD__(value='', type=str),
        }),
        ('bwa_QA', {
            'mem': __RECORD__(value=1250, type=int),
            'time': __RECORD__(value=235, type=int),
            'threads': __RECORD__(value=0, type=int),
            'conda': __RECORD__(value='', type=str),
        }),
        ('coverage_QA', {
            'mem': __RECORD__(value=1250, type=int),
            'time': __RECORD__(value=235, type=int),
            'conda': __RECORD__(value=f'{VPIPE_BASEDIR}/envs/smallgenomeutilities.yaml', type=str),
        }),
        ('msa', {
            'mem': __RECORD__(value=10000, type=int),
            'time': __RECORD__(value=235, type=int),
            'threads': __RECORD__(value=0, type=int),
            'conda': __RECORD__(value='', type=str),
        }),
        ('convert_to_ref', {
            'mem': __RECORD__(value=8000, type=int),
            'time': __RECORD__(value=235, type=int),
            'conda': __RECORD__(value=f'{VPIPE_BASEDIR}/envs/smallgenomeutilities.yaml', type=str),
        }),
        ('ref_bwa_index', {
            'mem': __RECORD__(value=2000, type=int),
            'time': __RECORD__(value=235, type=int),
            'conda': __RECORD__(value=f'{VPIPE_BASEDIR}/envs/bwa_align.yaml', type=str),
        }),
        ('bwa_align', {
            'mem': __RECORD__(value=1250, type=int),
            'time': __RECORD__(value=235, type=int),
            'threads': __RECORD__(value=0, type=int),
            'conda': __RECORD__(value='', type=str),

            'extra': __RECORD__(value='', type=str),
        }),
        ('ref_bowtie_index', {
            'mem': __RECORD__(value=2000, type=int),
            'time': __RECORD__(value=235, type=int),
            'conda': __RECORD__(value=f'{VPIPE_BASEDIR}/envs/bowtie_align.yaml', type=str),
        }),
        ('bowtie_align', {
            'mem': __RECORD__(value=1250, type=int),
            'time': __RECORD__(value=235, type=int),
            'threads': __RECORD__(value=0, type=int),
            'conda': __RECORD__(value='', type=str),

            'phred': __RECORD__(value='--phred33', type=str),
            'preset': __RECORD__(value='--local --sensitive-local', type=str),
            'extra': __RECORD__(value='', type=str),
        }),
        ('consensus_sequences', {
            'mem': __RECORD__(value=1250, type=int),
            'time': __RECORD__(value=235, type=int),
            'conda': __RECORD__(value=f'{VPIPE_BASEDIR}/envs/smallgenomeutilities.yaml', type=str),

            'min_coverage': __RECORD__(value=50, type=int),
            'qual_thrd': __RECORD__(value=15, type=int),
            'min_freq': __RECORD__(value=0.05, type=float),
        }),
        ('minor_variants', {
            'mem': __RECORD__(value=1000, type=int),
            'time': __RECORD__(value=235, type=int),
            'threads': __RECORD__(value=0, type=int),
            'conda': __RECORD__(value=f'{VPIPE_BASEDIR}/envs/smallgenomeutilities.yaml', type=str),

            'min_coverage': __RECORD__(value=100, type=int),
            'frequencies': __RECORD__(value=False, type=bool),
        }),
        ('coverage_intervals', {
            'mem': __RECORD__(value=1000, type=int),
            'time': __RECORD__(value=60, type=int),
            'threads': __RECORD__(value=0, type=int),
            'conda': __RECORD__(value=f'{VPIPE_BASEDIR}/envs/smallgenomeutilities.yaml', type=str),

            'overlap': __RECORD__(value=False, type=bool),
            'coverage': __RECORD__(value=50, type=int),
            'liberal': __RECORD__(value=True, type=bool),
        }),
        ('snv', {
            'mem': __RECORD__(value=10000, type=int),
            'time': __RECORD__(value=2880, type=int),
            'threads': __RECORD__(value=0, type=int),
            'conda': __RECORD__(value='', type=str),

            'alpha': __RECORD__(value=0.1, type=float),
            'ignore_indels': __RECORD__(value=False, type=bool),
            'coverage': __RECORD__(value=0, type=int),
            'shift': __RECORD__(value=3, type=int),
            'keep_files': __RECORD__(value=False, type=bool),
        }),
        ('lofreq', {
            'mem': __RECORD__(value=2000, type=int),
            'time': __RECORD__(value=60, type=int),
            'conda': __RECORD__(value='', type=str),

            'extra': __RECORD__(value='', type=str),
        }),
        ('aggregate', {
            'mem': __RECORD__(value=2000, type=int),
            'time': __RECORD__(value=235, type=int),
        }),
        ('haploclique', {
            'mem': __RECORD__(value=10000, type=int),
            'time': __RECORD__(value=1435, type=int),
            'conda': __RECORD__(value='', type=str),

            'relax': __RECORD__(value=True, type=bool),
            'no_singletons': __RECORD__(value=True, type=bool),
            'no_prob0': __RECORD__(value=True, type=bool),
            'clique_size_limit': __RECORD__(value=3, type=int),
            'max_num_cliques': __RECORD__(value=10000, type=int),
        }),
        ('haploclique_visualization', {
            'mem': __RECORD__(value=2000, type=int),
            'time': __RECORD__(value=235, type=int),
            'conda': __RECORD__(value=f'{VPIPE_BASEDIR}/envs/smallgenomeutilities.yaml', type=str),

            'region_start': __RECORD__(value=0, type=int),
            'region_end': __RECORD__(value=9719, type=int),
            'msa': __RECORD__(value='', type=str),
        }),
        ('savage', {
            'mem': __RECORD__(value=10000, type=int),
            'time': __RECORD__(value=1435, type=int),
            'threads': __RECORD__(value=0, type=int),
            'conda': __RECORD__(value='', type=str),

            'split': __RECORD__(value=20, type=int),
        }),
        ('web_visualization', {
            'mem': __RECORD__(value=2000, type=int),
            'time': __RECORD__(value=235, type=int),
            'conda': __RECORD__(value=f'{VPIPE_BASEDIR}/envs/visualization.yaml', type=str),
        }),
    ])

    def __init__(self):
        self.__members = {}

        vpipe_configfile = configparser.ConfigParser()
        vpipe_configfile.read('vpipe.config')

        for (section, properties) in VpipeConfig.__MEMBER_DEFAULT__.items():
            self.__members[section] = {}

            for (value, defaults) in properties.items():
                try:
                    if defaults.type == int:
                        cur_value = vpipe_configfile.getint(section, value)
                    elif defaults.type == float:
                        cur_value = vpipe_configfile.getfloat(section, value)
                    elif defaults.type == bool:
                        cur_value = vpipe_configfile.getboolean(section, value)
                    else:
                        cur_value = vpipe_configfile.get(section, value)
                    vpipe_configfile.remove_option(section, value)
                    state = 'user'
                except (configparser.NoSectionError, configparser.NoOptionError):
                    if value == 'threads' and section != 'general':
                        cur_value = defaults.value if defaults.value else self.__members[
                            'general']['threads']
                    elif value == 'conda':
                        cur_value = f'{VPIPE_BASEDIR}/envs/{section}.yaml' if len(
                            defaults.value) == 0 else defaults.value
                    else:
                        cur_value = defaults.value
                    state = 'DEFAULT'
                except ValueError as err:
                    raise ValueError("ERROR: Property '{}' of section '{}' has to be of type '{}', whereas you gave '{}'!".format(
                        value, section, defaults.type.__name__, vpipe_configfile[section][value])) from err

                if VPIPE_CONFIG_OPTS:
                    LOGGER.info(f"({state}) \t\t {section}: {value} = {cur_value}")

                self.__members[section][value] = cur_value

            if vpipe_configfile.has_section(section):
                if vpipe_configfile.items(section):
                    raise ValueError(
                        f"ERROR: Unrecognized options in section {section}: "
                        + ", ".join([option for option, _ in vpipe_configfile.items(section)]))
                vpipe_configfile.remove_section(section)

        sections_left = {section for section, _ in vpipe_configfile.items()} - {'DEFAULT'}
        if sections_left:
            raise ValueError(
                f"ERROR: Unrecognized sections in config file: "
                + ", ".join(sections_left))

    def __getattr__(self, name):
        try:
            return self.__members[name]
        except KeyError as err:
            raise ValueError(
                "ERROR: Section '{}' is not a valid section!".format(name)) from err


config = VpipeConfig()


# 2. glob patients/samples + store as TSV if file is not provided

# if file containing samples exists, proceed
# to build list of target files
if not os.path.isfile(config.input['samples_file']):
    # sample file does not exist, have to first glob
    # all patients' data and then construct sample
    # list that would pass QA checks

    patient_sample_pairs = glob_wildcards(
        "{}/{{patient_date}}/raw_data/{{file}}".format(config.input['datadir']))

    if not set(patient_sample_pairs.patient_date):
        LOGGER.warning(
            f"WARNING: No samples found in {config.input['datadir']}. Not generating {config.input['samples_file']}.")
    else:
        with open(config.input['samples_file'], 'w') as outfile:
            for i in set(patient_sample_pairs.patient_date):
                (patient, date) = [x.strip() for x in i.split("/") if x.strip()]
                outfile.write('{}\t{}\n'.format(patient, date))

    # TODO: have to preprocess patient files to filter likely failures
    # 1.) Determine 5%/95% length of FASTQ files
    # 2.) Determine whether FASTQ would survive


# 3. load patients from TSV and create list of samples
#
# This list is reused on subsequent runs

patient_list = []
patient_dict = {}
patient_record = typing.NamedTuple(
    "patient_record", [('patient_id', str), ('date', str)])

if not os.path.isfile(config.input['samples_file']):
    LOGGER.warning(
        f"WARNING: Sample list file {config.input['samples_file']} not found.")
else:
    with open(config.input['samples_file'], newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t')

        for row in spamreader:
            assert len(row) >= 2, "ERROR: Line '{}' does not contain at least two entries!".format(
                                  spamreader.line_num)
            patient_tuple = patient_record(patient_id=row[0], date=row[1])
            patient_list.append(patient_tuple)

            assert config.input['trim_percent_cutoff'] > 0 and config.input['trim_percent_cutoff'] < 1, "ERROR: 'trim_percent_cutoff' is expected to be a fraction (between 0 and 1), whereas 'trim_percent_cutoff'={}".format(
                config.input['trim_percent_cutoff'])
            assert patient_tuple not in patient_dict, "ERROR: sample '{}-{}' is not unique".format(
                row[0], row[1])

            if len(row) == 2:
                # All samples are assumed to have same read length and the default, 250 bp
                patient_dict[patient_tuple] = 250

            elif len(row) >= 3:
                # Extract read length from input.samples_file. Samples may have
                # different read lengths. Reads will be filtered out if read length
                # after trimming is less than trim_cutoff * read_length.
                patient_dict[patient_tuple] = int(row[2])

# 4. generate list of target files
all_files = []
alignments = []
vicuna_refs = []
references = []
consensus = []
trimmed_files = []
fastqc_files =[]
results = []
visualizations = []
datasets = []
IDs = []
for p in patient_list:

    alignments.append(
        "{sample_dir}/{patient}/{date}/alignments/REF_aln.bam".format(sample_dir=config.input['datadir'], patient=p.patient_id, date=p.date))
    if config.output['QA']:
        alignments.append(
            "{sample_dir}/{patient}/{date}/QA_alignments/coverage_ambig.tsv".format(sample_dir=config.input['datadir'], patient=p.patient_id, date=p.date))
        alignments.append(
            "{sample_dir}/{patient}/{date}/QA_alignments/coverage_majority.tsv".format(sample_dir=config.input['datadir'], patient=p.patient_id, date=p.date))

    vicuna_refs.append(
        "{sample_dir}/{patient}/{date}/references/vicuna_consensus.fasta".format(sample_dir=config.input['datadir'], patient=p.patient_id, date=p.date))
    references.append(
        "{sample_dir}/{patient}/{date}/references/ref_".format(sample_dir=config.input['datadir'], patient=p.patient_id, date=p.date))

    consensus.append(
        "{sample_dir}/{patient}/{date}/references/ref_ambig.fasta".format(sample_dir=config.input['datadir'], patient=p.patient_id, date=p.date))
    consensus.append(
        "{sample_dir}/{patient}/{date}/references/ref_majority.fasta".format(sample_dir=config.input['datadir'], patient=p.patient_id, date=p.date))

    trimmed_files.append(
        "{sample_dir}/{patient}/{date}/preprocessed_data/R1.fastq.gz".format(sample_dir=config.input['datadir'], patient=p.patient_id, date=p.date))
    if config.input['paired']:
        trimmed_files.append("{sample_dir}/{patient}/{date}/preprocessed_data/R2.fastq.gz".format(
            sample_dir=config.input['datadir'], patient=p.patient_id, date=p.date))

    fastqc_files.append(
        "{sample_dir}/{patient}/{date}/extracted_data/R1_fastqc.html".format(sample_dir=config.input['datadir'], patient=p.patient_id, date=p.date))
    if config.input['paired']:
        fastqc_files.append("{sample_dir}/{patient}/{date}/extracted_data/R2_fastqc.html".format(
            sample_dir=config.input['datadir'], patient=p.patient_id, date=p.date))

    datasets.append("{sample_dir}/{patient}/{date}".format(
        sample_dir=config.input['datadir'], patient=p.patient_id, date=p.date))
    IDs.append(('{}-{}').format(p.patient_id, p.date))

    # SNV
    if config.output['snv']:
        if config.general['snv_caller'] == 'shorah':
            results.append("{sample_dir}/{patient}/{date}/variants/SNVs/snvs.csv".format(
                sample_dir=config.input['datadir'], patient=p.patient_id, date=p.date))
        elif config.general['snv_caller'] == 'lofreq':
            results.append("{sample_dir}/{patient}/{date}/variants/SNVs/snvs.vcf".format(
                sample_dir=config.input['datadir'], patient=p.patient_id, date=p.date))
    # local haplotypes
    if config.output['local']:
        results.append(
            "{sample_dir}/{patient}/{date}/variants/SNVs/snvs.csv".format(sample_dir=config.input['datadir'], patient=p.patient_id, date=p.date))
    # global haplotypes
    if config.output['global']:
        if config.general['haplotype_reconstruction'] == 'savage':
            results.append("{sample_dir}/{patient}/{date}/variants/global/contigs_stage_c.fasta".format(
                sample_dir=config.input['datadir'], patient=p.patient_id, date=p.date))
        elif config.general['haplotype_reconstruction'] == 'haploclique':
            results.append("{sample_dir}/{patient}/{date}/variants/global/quasispecies.bam".format(
                sample_dir=config.input['datadir'], patient=p.patient_id, date=p.date))

    # visualization
    if config.output['visualization']:
        visualizations.append("{sample_dir}/{patient}/{date}/visualization/index.html".format(
            sample_dir=config.input['datadir'], patient=p.patient_id, date=p.date))

    # merge lists contaiing expected output
    all_files = alignments + consensus + results + visualizations

IDs = ','.join(IDs)

# 5. Locate reference and parse reference identifier
if not VPIPE_BENCH:
    reference_file = config.input['reference']
    if not os.path.isfile(reference_file):
        reference_file_alt = os.path.join("references", reference_file)
        LOGGER.warning(
            f"WARNING: Reference file {reference_file} not found. Trying {reference_file_alt}.")
        reference_file = reference_file_alt
        if not os.path.isfile(reference_file):
            raise ValueError(
                f"ERROR: Reference file {reference_file} not found.")

    with open(reference_file, 'r') as infile:
        reference_name = infile.readline().rstrip()
    reference_name = reference_name.split('>')[1]
    reference_name = reference_name.split(' ')[0]

else:
    try:
        reference_name
    except NameError:
        reference_name = ''
