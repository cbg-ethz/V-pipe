import csv
import collections
import configparser
from typing import Dict, Any, NamedTuple
import os

__author__ = "Susana Posada-Cespedes"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"

# 1. Parse config file

# Import VpipeConfig class defining defaults
include: "config_default.smk"


def _deep_merge(
        dict1: Dict[str, Dict[str, Any]], dict2: Dict[str, Dict[str, Any]]
) -> Dict[str, Dict[str, Any]]:
    dict_out = dict1.copy()
    for key, value in dict2.items():
        dict_out[key] = {**dict1.get(key, {}), **value}
    return dict_out


class VpipeBenchConfig(VpipeConfig):
    'Class used to encapsulate the configuration properties used by benchmarking scripts'

    __RECORD__ = VpipeConfig.__RECORD__

    __MEMBER_DEFAULT__ = _deep_merge(
        VpipeConfig.__MEMBER_DEFAULT__,
        collections.OrderedDict([
            ('general', {
                'seed': __RECORD__(value=42, type=int),
                'simulate': __RECORD__(value=True, type=bool),
            }),
            ('applications', {
                'simBench': __RECORD__(
                    value=f"{VPIPE_BASEDIR}/scripts/simBench.py", type=str),
                'art': __RECORD__(value="art_illumina", type=str),
                'testBench': __RECORD__(
                    value=f"{VPIPE_BASEDIR}/scripts/testBench.py", type=str),
                'coverageIntervals' : __RECORD__(
                    value=f"{VPIPE_BASEDIR}/scripts/coverageIntervals.py", type=str),
                'snakemake': __RECORD__(value="snakemake", type=str),
            }),
            ('benchmark', {
                'aligners': __RECORD__(value='ngshmmalign,bwa,bowtie', type=str),
                'snv_callers': __RECORD__(value='shorah,lofreq', type=str),
                'snakemake_options': __RECORD__(value='--cores 4 -p', type=str)
            }),
            ('simulate_master', {
                'mem': __RECORD__(value=2000, type=int),
                'time': __RECORD__(value=30, type=int),
                'conda': __RECORD__(
                    value=f'{VPIPE_BASEDIR}/envs/simbench.yaml', type=str),

                'genome_length': __RECORD__(value=3000, type=int),
                'seed': __RECORD__(value=0, type=int),
            }),
            ('simulate_haplotypes', {
                'mem': __RECORD__(value=2000, type=int),
                'time': __RECORD__(value=30, type=int),
                'conda': __RECORD__(
                    value=f'{VPIPE_BASEDIR}/envs/simbench.yaml', type=str),
            }),
            ('simulate_reads', {
                'mem': __RECORD__(value=2000, type=int),
                'time': __RECORD__(value=30, type=int),
                'conda': __RECORD__(
                    value=f'{VPIPE_BASEDIR}/envs/simbench.yaml', type=str),

                'num_reads': __RECORD__(value=False, type=bool),
                'high_quality': __RECORD__(value=True, type=bool),
            }),
            ('run_simBench', {
                'mem': __RECORD__(value=5000, type=int),
                'time': __RECORD__(value=90, type=int),
            }),
            ('test_snv', {
                'mem': __RECORD__(value=2000, type=int),
                'time': __RECORD__(value=60, type=int),
                'conda': __RECORD__(value='', type=str),

                're_msa': __RECORD__(value=False, type=bool),
            }),
            ('aggregate', {
                'mem': __RECORD__(value=2000, type=int),
                'time': __RECORD__(value=235, type=int),
            }),
            ('run_vpipeBench', {
                'mem': __RECORD__(value=5000, type=int),
                'time': __RECORD__(value=1440, type=int),
            }),
        ])
    )


VPIPE_CONFIG = VpipeBenchConfig

include: "common.smk"


# 2. Parse file containing info about simulated data sets
sample_dict = {}
sample_record = NamedTuple(
    "sample_record", [('sample_name', str), ('date', str)])
datasets = []

if not os.path.isfile(config.input['samples_file']):
    raise ValueError(
        f"ERROR: Sample list file {config.input['samples_file']} not found.")
else:
    with open(config.input['samples_file'], newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t')

        for row in spamreader:
            assert len(row) <= 16, "ERROR: Line '{}' contains more entries than expected".format(
                spamreader.line_num)
            sample_tuple = sample_record(sample_name=row[0], date=row[1])

            datasets.append("{sample_dir}/{ID}/{date}".format(
                sample_dir=config.input['datadir'], ID=row[0], date=row[1]))

            if len(row) == 2:
                # All data sets are assumed to have default specifications
                sample_dict[sample_tuple] = {"read_len": 250, "haplotype_seqs": None, "num_haplotypes": 5, "coverage": 500, "fragment_size": '600,100',
                                             "freq_dstr": 'geom', "freq_param": 0.75, "mut_rate": 0.1, "del_rate": 0.0, "ins_rate": 0.0, "no_frameshifts": True, "del_len": ''}

            else:
                # All other features (except seed) must be specified
                sample_dict[sample_tuple] = {"read_len": int(row[2]), "haplotype_seqs": row[3], "num_haplotypes": int(row[4]), "coverage": int(
                    row[5]), "fragment_size": row[6], "freq_dstr": row[7], "freq_param": row[8], "mut_rate": float(row[9]), "del_rate": float(row[10]), "ins_rate": float(row[11]), "no_frameshifts": row[12], "del_len": row[13]}

            if len(row) == 15:
                # Seed is specified
                sample_dict[sample_tuple]['seed'] = int(row[14])
            else:
                # Seed is not provided / note that potential replicates would have the same seed
                sample_dict[sample_tuple]['seed'] = config.general['seed']


# 3. V-pipe expects a reference as input. We need to "mask" this behaviour
if config.input['reference']:
    # Locate reference file
    if os.path.isfile(config.input['reference']):
        reference_file = config.input['reference']
        reference_name = get_reference_name(reference_file)
    elif os.path.isfile(os.path.join("references", config.input['reference'])):
        reference_file = os.path.join(
            "references", config.input['reference'])
        reference_name = get_reference_name(reference_file)
    else:
        # If reference file not found, create it
        reference_file = "references/haplotype_master.fasta"
        reference_name = "master"
else:
    reference_file = "references/haplotype_master.fasta"
    reference_name = "master"


# Auxiliary functions


def get_haplotype_seqs(wildcards):
    sample_tuple = sample_record(
        sample_name=wildcards.sample_name, date=wildcards.date)
    haplotype_seqs = sample_dict[sample_tuple]['haplotype_seqs']
    val = ''
    if haplotype_seqs is not None:
        if len(haplotype_seqs) > 0 and haplotype_seqs.upper() not in ["NA", "N/A"]:
            val = haplotype_seqs
    return val


def get_no_FR(wildcards):
    sample_tuple = sample_record(
        sample_name=wildcards.sample_name, date=wildcards.date)
    no_FR = sample_dict[sample_tuple]['no_frameshifts']
    return '-fr' if no_FR else ''


def get_del_len(wildcards):
    sample_tuple = sample_record(
        sample_name=wildcards.sample_name, date=wildcards.date)
    del_len = sample_dict[sample_tuple]['del_len']
    del_len = del_len.strip()
    if not del_len or del_len.upper() in ["NA", "N/A"]:
        val = ''
    else:
        val = '-dl {}'.format(del_len)
    return val


def get_freq_params(wildcards):
    sample_tuple = sample_record(
        sample_name=wildcards.sample_name, date=wildcards.date)
    freq_dstr = sample_dict[sample_tuple]['freq_dstr']
    freq_param = sample_dict[sample_tuple]['freq_param']
    if freq_dstr == 'geom':
        val = '-gr {}'.format(freq_param)
    elif freq_dstr == 'dirichlet':
        val = '-dc {}'.format(freq_param)
    else:
        val = ''
    return val


def get_freq_aux(wildcards):
    sample_tuple = sample_record(
        sample_name=wildcards.sample_name, date=wildcards.date)
    freq_dstr = sample_dict[sample_tuple]['freq_dstr']
    freq_param = sample_dict[sample_tuple]['freq_param']
    if freq_dstr == 'geom':
        val = '-gr {}'.format(freq_param)
    elif freq_dstr == 'dirichlet':
        infile = os.path.join(wildcards.sample_dir, wildcards.sample_name,
                              wildcards.date, "references/haplotypes/haplotype_frequencies.fasta")
        val = '-df {}'.format(infile)
    else:
        val = ''
    return val


def input_snv(wildcards):
    input = []
    if config.general['snv_caller'] == 'shorah':
        input.append(os.path.join(wildcards.sample_dir, wildcards.sample_name,
                                  wildcards.date, "variants/SNVs/snvs.csv"))
        input.append(os.path.join(wildcards.sample_dir, wildcards.sample_name,
                                  wildcards.date, "variants/coverage_intervals.tsv"))
    elif config.general['snv_caller'] == 'lofreq':
        input.append(os.path.join(wildcards.sample_dir, wildcards.sample_name,
                                  wildcards.date, "variants/SNVs/snvs.vcf"))
        input.append(os.path.join("variants/coverage.tsv"))
    return input
