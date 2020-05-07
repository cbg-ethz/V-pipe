import collections
import configparser
import os
import typing

__author__ = "Susana Posada-Cespedes"
__author__ = "David Seifert"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"


# Class to parse config file
#
# the config object sets defaults for a number of parameters and
# validates parameters given by the user to be sensible

VPIPE_CONFIG_OPTS = True if os.environ.get(
    'VPIPE_CONFIG_OPTS') is not None else False

class _SectionWrapper(object):
    def __init__(self, config: 'VpipeConfig', section):
        if not config.has_section(section):
            raise KeyError(
                "ERROR: Section '{}' is not a valid section!".format(section))
        self._config = config
        self._section = section

    def __setitem__(self, option, value):
        self._config.set_option(self._section, option, value)

    def __getitem__(self, option):
        return self._config.get_option(self._section, option)


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
        }),
        ('output', {
            'QA': __RECORD__(value=False, type=bool),
            'snv': __RECORD__(value=True, type=bool),
            'local': __RECORD__(value=True, type=bool),
            'global': __RECORD__(value=True, type=bool),
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

            'ref_panel': __RECORD__(value='references/5-Virus-Mix.fasta', type=str),
        }),
        ('coverage_QA', {
            'mem': __RECORD__(value=1250, type=int),
            'time': __RECORD__(value=235, type=int),
            'conda': __RECORD__(value=f'{VPIPE_BASEDIR}/envs/smallgenomeutilities.yaml', type=str),

            'target': __RECORD__(value='HXB2:6614-6812,7109-7217,7376-7478,7601-7634', type=str),
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
            'maxins': __RECORD__(value=None, type=int),
            'extra': __RECORD__(value='', type=str),
        }),
        ('consensus_sequences', {
            'mem': __RECORD__(value=1250, type=int),
            'time': __RECORD__(value=235, type=int),
            'conda': __RECORD__(value=f'{VPIPE_BASEDIR}/envs/smallgenomeutilities.yaml', type=str),

            'min_coverage': __RECORD__(value=50, type=int),
            'n_coverage': __RECORD__(value=5, type=int),
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

            'consensus': __RECORD__(value=True, type=bool),
            'alpha': __RECORD__(value=0.1, type=float),
            'ignore_indels': __RECORD__(value=False, type=bool),
            'coverage': __RECORD__(value=0, type=int),
            'shift': __RECORD__(value=3, type=int),
            'keep_files': __RECORD__(value=False, type=bool),
        }),
        ('samtools_index', {
            'mem': __RECORD__(value=2000, type=int),
            'time': __RECORD__(value=20, type=int),
            'conda': __RECORD__(value=f'{VPIPE_BASEDIR}/envs/lofreq.yaml', type=str),
        }),
        ('lofreq', {
            'mem': __RECORD__(value=2000, type=int),
            'time': __RECORD__(value=60, type=int),
            'conda': __RECORD__(value='', type=str),

            'consensus': __RECORD__(value=True, type=bool),
            'extra': __RECORD__(value='', type=str),
        }),
        ('alignment_coverage', {
            'mem': __RECORD__(value=1000, type=int),
            'time': __RECORD__(value=60, type=int),
            'conda': __RECORD__(value=f'{VPIPE_BASEDIR}/envs/smallgenomeutilities.yaml', type=str),

            'coverage': __RECORD__(value=5, type=int),
        }),
        ('stats', {
            'conda': __RECORD__(value=f'{VPIPE_BASEDIR}/envs/sam2bam.yaml', type=str),
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
        })
    ])

    def __init__(self):
        # track all options (explicitly-set and defaults)
        self.__members = {}

        # track exclicitly-set options
        self._vpipe_configfile = configparser.ConfigParser()
        self._vpipe_configfile.read('vpipe.config')

        # clone to validate vpipe.config file 
        vpipe_configfile = configparser.ConfigParser()
        vpipe_configfile.read('vpipe.config')

        for (section, properties) in self.__MEMBER_DEFAULT__.items():
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
                    LOGGER.info(
                        f"({state}) \t\t {section}: {value} = {cur_value}")

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

    def __getattr__(self, section_name) -> _SectionWrapper:
        return _SectionWrapper(self, section_name)

    def has_section(self, section):
        return section in self.__members

    def set_option(self, section, option, value):
        if section not in self.__members:
            raise KeyError(
                "ERROR: Section '{}' is not a valid section!".format(section))
        self.__members[section][option] = value

        if section not in self._vpipe_configfile:
            self._vpipe_configfile[section] = {}
        # now add this explicitly set option
        self._vpipe_configfile[section][option] = value

    def get_option(self, section, option):
        if section not in self.__members:
            raise KeyError(
                "ERROR: Section '{}' is not a valid section!".format(section))
        elif option not in self.__members[section]:
            raise ValueError(
                "ERROR: Section '{}' has no property '{}'!".format(section, option))
        else:
            return self.__members[section][option]

    def write(self, outfile):
        self._vpipe_configfile.write(outfile)
