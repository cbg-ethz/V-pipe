#!/usr/bin/env python3

import os
import argparse
from alignmentIntervals import read_fasta

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sh
import numpy as np
import pandas as pd


__author__ = "Susana Posada-Cespedes"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"


DBG = True if os.environ.get('DBG') is not None else False


def parse_args():
    """ Set up the parsing of command-line arguments """

    parser = argparse.ArgumentParser(
        description="Benchmark: test",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        "-f", required=True, default=None, metavar='FASTA',
        dest='haplotype_seqs',
        help="Fasta file containing either the sequences of the true "
             "haplotypes or haplotypes sequences (msa) already reported using "
             "the same indexing as the reference/consensus sequence"
    )
    requiredNamed.add_argument(
        "-s", required=True, default=None, metavar='CSV', dest='snvs',
        help="File containing called SNVs"
    )
    requiredNamed.add_argument(
        "-N", required=False, default='sample', metavar='STR',
        dest='sampleID', help="Patient/sample identifiers"
    )
    parser.add_argument(
        "-m", required=False, default=None, metavar='FASTA',
        dest='haplotype_master', type=str,
        help="Fasta file containing the sequence with respect to which SNVs "
             "were called"
    )
    parser.add_argument(
        "--ref", required=False, default=None, metavar='FASTA',
        dest='reference', type=str,
        help="Fasta file containing the reference sequence with respect to "
             "which reads were aligned"
    )
    parser.add_argument(
        "-d", required=False, default='unif', metavar='str', dest='freq_dstr',
        type=str, choices=['unif', 'geom', 'dirichlet', 'cust'],
        help="Distribution of haplotype frequencies"
    )
    parser.add_argument(
        "-gr", required=False, default=0.75, metavar='FLOAT', dest='ratio',
        type=float, help="Sucess probability for the geometric distribution"
    )
    parser.add_argument(
        "-df", required=False, default=None, metavar='FASTA',
        dest='dirichlet_freqs', type=str,
        help="File containing haplotype frequencies"
    )
    parser.add_argument(
        "-ci", required=False, default=None, metavar='chrm:start-end',
        dest='coverage_intervals', type=str,
        help="File containing coverage intervals"
    )
    parser.add_argument(
        "--no-expansion", required=False, default=False, action='store_true',
        dest='no_expansion',
        help="Coverage intervals do not correspond to region use to run "
             "ShoRAH, but the actual target region"
    )
    parser.add_argument(
        "--caller", required=False, default='shorah', metavar='str',
        dest='snv_caller', type=str, choices=['shorah', 'lofreq'],
        help="Inidcate if other software different from ShoRAH was used for "
             "SNV calling"
    )
    parser.add_argument(
        "-wl", required=False, default=201, metavar='INT', dest='window_len',
        type=int,
        help="Window length used by ShoRAH to construct overlapping windows"
    )
    parser.add_argument(
        "-ws", required=False, default=3, metavar='INT', dest='window_shift',
        type=int,
        help="Number of window shifts used by ShoRAH to construct overlapping "
             "windows"
    )
    parser.add_argument(
        "-cf", required=False, default=None, metavar='TXT', dest='coverage',
        type=str,
        help="File to read coverage per window used by ShoRAH, or a "
             "tab-separated values file containing coverage per locus"
    )
    parser.add_argument(
        "-ms", required=False, default=False, action='store_true', dest='msa',
        help="Indicate if the multiple sequence alignment including "
             "reference/consensus sequence should be constructed"
    )
    parser.add_argument(
        "--only-dels", required=False, default=False, action='store_true',
        dest='only_deletions',
        help="Indicate if only performance based on deletions should reported"
    )
    parser.add_argument(
        "--long-dels", required=False, default=False, action='store_true',
        dest='long_deletions',
        help="Indicate if deletions should be parsed as multipe-base deletions"
    )
    parser.add_argument(
        "-t", required=False, default=False, action='store_true',
        dest='output_true',
        help="Indicate if file containing expected SNVs should be reported. "
             "Report using 1-based indexing for the position"
    )
    parser.add_argument(
        "-mafft", required=False, default="mafft", metavar='PATH',
        dest='mafft', type=str,
        help="Path to binaries for the multiple sequence aligner MAFFT"
    )
    parser.add_argument(
        "-of", required=False, default='performance.tsv', metavar='OUTPUT',
        dest='outfile', type=str,
        help="Output file - file containing expected SNVs"
    )
    parser.add_argument(
        "-od", required=False, default=None, metavar='DIR', dest='outdir',
        type=str, help="Output directory for intermediate files"
    )

    return parser.parse_args()


def frequencies(freq_dstr, num_haplotypes, ratio=0.75, infile=None):
    "Compute the expected haplotype frequencies"
    if freq_dstr == 'unif':
        haplotype_freqs = np.repeat(1 / num_haplotypes, num_haplotypes)
    elif freq_dstr == 'geom':
        haplotype_freqs = [ratio**(i + 1) for i in range(num_haplotypes)]
        haplotype_freqs = np.asarray(haplotype_freqs)
        haplotype_freqs = haplotype_freqs / np.sum(haplotype_freqs)
    elif freq_dstr == 'dirichlet':
        # Read haplotype frequencies from output file
        if infile is None:
            raise IOError(
                    "Input file containing haplotype frequencies is expected")
        ids, haplotype_freqs = read_fasta(infile)
        haplotype_freqs = np.asarray(haplotype_freqs, dtype=float)

    return haplotype_freqs


def parse_info(df, snvcaller):
    if snvcaller == 'shorah':
        df_info = pd.DataFrame.from_dict(
            [dict([entry.strip().split("=") for entry in line.split(";")])
                for line in df["INFO"]]).astype('float')
        # We ignore columns with 0-counts to compute the SNV frequency. A zero
        # count means that the SNV was not found in the corresponding window.
        df_freq = df_info[["Freq1", "Freq2", "Freq3"]].copy()
        df_freq[df_freq == 0] = np.nan
        df_freq = df_freq.mean(axis=1)
    elif snvcaller == 'lofreq':
        df["INFO"] = df["INFO"].str.replace("INDEL", "INDEL=1")
        df_info = pd.DataFrame.from_dict(
            [dict([entry.strip().split("=") for entry in line.split(";")])
                for line in df["INFO"]])
        df_freq = df_info["AF"]
    return df_freq


def parse_vcf(snvfile, snvcaller):
    # Read VCF file to infer how many lines to skip
    skiplines = 0
    with open(snvfile, 'r') as infile:
        for line in infile:
            if not line.startswith('##'):
                break
            skiplines += 1
    try:
        df_snvs = pd.read_csv(snvfile, sep="\t", skiprows=skiplines, header=0,
                              compression=None)
        df_snvs = df_snvs.rename(columns={'#CHROM': 'CHROM'})
        df_snvs['FREQ'] = parse_info(df_snvs, snvcaller)
    except pd.errors.EmptyDataError:
        df_snvs = pd.DataFrame()
    return df_snvs


def true_snvs(haplotype_master_arr, haplotype_master, haplotype_seqs,
              num_haplotypes, haplotype_freqs, long_deletions, alphabet):
    """
    Extract expected SNVs using the MSA of the true haplotype sequences and
    the reference sequence
    """
    # loci = np.arange(haplotype_master_arr.size)
    haplotype_idx = np.arange(num_haplotypes)
    variants = haplotype_master_arr != haplotype_seqs

    df_snvs = pd.DataFrame(columns=('POS', 'REF', 'ALT', 'FREQ', 'HAPLOTYPES'))
    num_snvs = 0
    for locus in range(haplotype_master_arr.size):
        idxs = variants[:, locus]
        if np.any(idxs):
            var = haplotype_seqs[idxs, locus]
            snv_freq = haplotype_freqs[idxs]
            if np.sum(idxs) == 1:
                df_snvs.loc[num_snvs] = [
                    locus, haplotype_master_arr[locus].decode(),
                    var[0].decode(), snv_freq[0],
                    haplotype_idx[idxs].astype(str)[0]]
                num_snvs += 1
            else:
                for base in alphabet:
                    idxs_base = var == base
                    if np.sum(idxs_base) > 0:
                        hap_aux = ','.join(
                            haplotype_idx[idxs][idxs_base].astype(str))
                        df_snvs.loc[num_snvs] = [
                            locus,
                            haplotype_master_arr[locus].decode(),
                            base.decode(), np.sum(snv_freq[idxs_base]),
                            hap_aux]
                        num_snvs += 1
    df_snvs["POS"] = df_snvs["POS"].astype(int)

    if long_deletions:
        df_long_dels = pd.DataFrame({
            'POS': pd.Series([], dtype='int'),
            'REF': pd.Series([], dtype='str'),
            'ALT': pd.Series([], dtype='str'),
            'FREQ': pd.Series([], dtype='float'),
            'HAPLOTYPES': pd.Series([], dtype='str')})
        for idx, seq in enumerate(haplotype_seqs):
            is_deletion = np.concatenate(([0], seq == b'-', [0]))
            intervals = np.where(
                np.abs(np.diff(is_deletion)) == 1)[0].reshape(-1, 2)
            if intervals.size > 0:
                assert (intervals[:, 0] > 0).all(), (
                    "Deletion reported in the first reference position")
                # Deletions are by convention reported at the preceding
                # position
                dict_dels = {
                    'POS': intervals[:, 0] - 1,
                    'REF': [
                        haplotype_master[(x[0] - 1):x[1]] for x in intervals],
                    'ALT': [haplotype_master[x[0] - 1] for x in intervals],
                    'FREQ': [haplotype_freqs[idx]] * intervals.shape[0],
                    'HAPLOTYPES': [
                        str(haplotype_idx[idx])] * intervals.shape[0]
                }
                df_tmp = pd.DataFrame.from_dict(dict_dels)
                df_long_dels = pd.concat(
                    [df_long_dels, df_tmp], ignore_index=True)
        # Merge deletions found in different haplotypes together
        grpby = df_long_dels.set_index(["POS", "REF", "ALT"])[
            ["FREQ", "HAPLOTYPES"]].groupby(["POS", "REF", "ALT"])

        df_long_dels = pd.concat(
            [grpby["FREQ"].sum(),
             grpby["HAPLOTYPES"].apply(lambda s: ",".join(s))], axis=1)
        df_long_dels.reset_index(inplace=True)

        # Drop one-base deletions
        del_mask = df_snvs["ALT"].str.startswith('-')
        df_snvs = df_snvs[~del_mask]
        df_snvs = pd.concat(
                [df_snvs, df_long_dels], ignore_index=True)
        df_snvs = df_snvs.set_index(["POS", "REF", "ALT"])
        df_snvs = df_snvs.sort_index()
        df_snvs.reset_index(inplace=True)

    return df_snvs


def mafft(infile, outfile, max_iter=1000, thrd=4, mafft='mafft'):
    "Use MAFFT to obtain the multiple sequence alignment"
    # --nuc sequences are nucleotide
    # --localpair pairwise alignments
    # --maxiterate number of iterative refinement
    cmd = sh.Command(mafft)
    cmd = cmd.bake('--nuc')
    cmd = cmd.bake('--preservecase')
    cmd = cmd.bake('--maxiterate', max_iter)
    cmd = cmd.bake('--localpair')
    cmd = cmd.bake('--thread', thrd)
    cmd = cmd.bake(infile)
    cmd = cmd.bake(_out=outfile)

    print(cmd)
    cmd()


def consecutive(array, stepsize=1):
    return np.split(array, np.where(np.diff(array) != stepsize)[0] + 1)


def target_snvs(start_region, end_region, start_locus, long_deletions,
                end_locus=None):
    if long_deletions:
        is_contained = (start_locus >= start_region) & \
            (end_locus < end_region)
    else:
        is_contained = (start_locus >= start_region) & \
            (start_locus < end_region)
    return is_contained


def main():

    args = parse_args()

    alphabet = ['-', 'A', 'C', 'G', 'T']
    alphabet = np.array(alphabet, dtype='c')

    # Compute average frequency for SNVs called using ShoRAH
    df_snvs = parse_vcf(args.snvs, args.snv_caller)

    if df_snvs.empty:
        print("No called SNVs")
        with open(args.outfile, 'w') as outfile:
            outfile.write('ID\tTP\tFP\tFN\tTN\n')
        return

    # Drop insertions
    ins_mask = df_snvs["ALT"].str.len() > 1
    df_snvs = df_snvs[~ins_mask]

    if args.only_deletions:
        # Only look at deletions
        is_deletion = df_snvs["REF"].str.len() > 1
        df_snvs = df_snvs[is_deletion]

    if df_snvs.empty:
        print("No called SNVs")
        with open(args.outfile, 'w') as outfile:
            outfile.write('ID\tTP\tFP\tFN\tTN\n')
        return

    if not args.long_deletions:
        # Unroll deletions into one-base deletions
        del_mask = df_snvs["REF"].str.len() > 1
        assert (df_snvs.loc[del_mask, "ALT"] == df_snvs.loc[
            del_mask, "REF"].str[0]).all(), (
                "Reference base preceding deletion does not match")

        del_len = df_snvs.loc[del_mask, "REF"].str.len() - 1
        df_del = pd.DataFrame(
            np.repeat(df_snvs[del_mask].values, del_len.to_list(), axis=0))
        df_del.columns = df_snvs.columns
        df_del["ALT"] = '-'
        aux_idx = 0
        aux_pos = df_del.columns.get_loc("POS")
        aux_ref = df_del.columns.get_loc("REF")
        for idx, row in df_snvs[del_mask].iterrows():
            # ignore first base as it corresponds to the reference at the
            # preceding locus
            ref = list(row["REF"][1:])
            pos = [row["POS"] + x + 1 for x in range(len(ref))]
            df_del.iloc[aux_idx:(aux_idx + del_len[idx]), aux_pos] = pos
            df_del.iloc[aux_idx:(aux_idx + del_len[idx]), aux_ref] = ref
            aux_idx += del_len[idx]

        # Handle special case: reference sequence might contain a gap character
        # and a long deletion could include it. When unrolling long deletions
        # the REF and ALT fields will contain both gaps symbols
        is_gap = (df_del["REF"] == '-') & (df_del["ALT"] == '-')
        df_del = df_del[~is_gap]

        # Remove previous rows corresponding to deletions and add the one-base
        # deletions
        df_snvs = df_snvs[~del_mask]
        df_snvs = pd.concat(
                [df_snvs, df_del], ignore_index=True)
        df_snvs = df_snvs.set_index(["POS", "ALT", "REF"])
        df_snvs = df_snvs.sort_index()

        # Merge on POS and ALT
        grpby = df_snvs.set_index("CHROM", append=True)[
            ["INFO", "FREQ"]].groupby(["POS", "ALT", "REF", "CHROM"])
        df_snvs = pd.concat([grpby["INFO"].apply(lambda s: ";".join(s)),
                             grpby["FREQ"].sum()], axis=1)
        # grpby["REF"].first() # If not part of the index

    outdir = args.outdir if args.outdir is not None else os.getcwd()
    if args.haplotype_master is not None:
        # Parse file containing reference/consensus sequence (sequence w.r.t
        # which SNVs were called)
        header, haplotype_master = read_fasta(args.haplotype_master)
        header = header[0]
        haplotype_master = haplotype_master[0].upper()
        haplotype_master_array = np.array(list(haplotype_master))
        reference_len = haplotype_master_array.size

        if args.msa:
            # Expected if cohort consensus has gaps
            if args.reference:
                tmp, reference = read_fasta(args.reference)
                reference = reference[0].upper()
                reference = np.array(list(reference))
                assert reference.size == haplotype_master_array.size, (
                    "Reference and cohort consensus have different lengths")
                idxs_gaps = haplotype_master_array == '-'
                haplotype_master_array[idxs_gaps] = reference[idxs_gaps]
                args.haplotype_master = os.path.join(outdir,
                                                     'cohort_consensus.fasta')
                cohort_consensus = SeqRecord(Seq(''.join(
                    haplotype_master_array)), id=header, description="")
                with open(args.haplotype_master, 'w') as outfile:
                    SeqIO.write(cohort_consensus, outfile, "fasta")

            haplotype_master_array = haplotype_master_array.astype('c')
            # construct msa: haplotypes + reference/consensus sequence
            infile = os.path.join(outdir, "tmp.fasta")
            sh.cat([args.haplotype_seqs, args.haplotype_master], _out=infile)
            msa_file = os.path.join(outdir, 'haplotypes_re-msa.fasta')
            mafft(infile, msa_file, mafft=args.mafft)
            os.remove(infile)
            # Parse fasta file containing msa
            haplotype_ids, haplotype_seqs = read_fasta(msa_file)
            num_haplotypes = len(haplotype_ids) - 1

            haplotype_ref = haplotype_seqs[-1]
            haplotype_ref = haplotype_ref.upper()
            haplotype_ref = np.array(haplotype_ref, dtype='c')
            if haplotype_ref.size != reference_len:
                assert haplotype_ref.size > reference_len, (
                    "Length of the consensus/reference sequence after the "
                    "MSA is smaller")
                # Deletions '-' were placed on the consensus/reference
                # sequence after the msa
                idx_master = 0
                idx_ref = 0
                idxs_ref = np.arange(haplotype_ref.size)
                del_idxs = np.zeros(haplotype_ref.size, dtype=bool)
                for i in range(haplotype_ref.size - reference_len):
                    left = min(reference_len + i - idx_ref,
                               haplotype_master_array[idx_master:].size)
                    idxs = haplotype_ref[idx_ref:(
                        idx_ref + left)] == haplotype_master_array[idx_master:]
                    aux = idxs_ref[idx_ref:(idx_ref + left)][~idxs]
                    if aux.size == 0:
                        # gaps '-' were placed until the end of haplotype_ref
                        del_idxs[(idx_ref + left):] = True
                        break
                    else:
                        idx_master = aux[0] - i
                        idx_ref = aux[0] + 1
                        del_idxs[aux[0]] = True

                assert np.all(
                    haplotype_ref[~del_idxs] == haplotype_master_array
                    ), "After substracting gaps sequences do not agree"
                assert np.all(
                    haplotype_ref[del_idxs] == b'-'
                    ), "All substracted loci do not correspond to '-'"

            # Parse sequences of the true haplotype
            haplotype_ids = haplotype_ids[0:num_haplotypes]
            haplotype_seqs = haplotype_seqs[0:num_haplotypes]
            haplotype_seqs_array = np.array(haplotype_seqs, dtype='c')
            # Remove insertions with respect to consensus/reference sequence
            if haplotype_ref.size != reference_len:
                haplotype_seqs_array = haplotype_seqs_array[:, ~del_idxs]
            # Restore gaps into the master sequence
            if args.reference:
                haplotype_master_array[idxs_gaps] = b'-'
        else:
            # Sequences of true haplotypes are already reported using the same
            # indexing as reference/consensus
            # Parse file containing true haplotype sequences
            haplotype_ids, haplotype_seqs = read_fasta(args.haplotype_seqs)
            num_haplotypes = len(haplotype_ids)
            haplotype_seqs_array = np.array(haplotype_seqs, dtype='c')
            haplotype_master_array = haplotype_master_array.astype('c')
    else:
        # if master sequence is not provided, report with respect to the
        # consensus. Note that SNVs are called with respect to the cohort
        # consensus.
        from scipy.stats import mode
        outfile = os.path.join(outdir, 'true_haplotype_msa.fasta')
        mafft(args.haplotype_seqs, outfile, mafft=args.mafft)
        haplotype_ids, haplotype_seqs = read_fasta(outfile)
        num_haplotypes = len(haplotype_ids)
        haplotype_seqs_array = np.array(haplotype_seqs, dtype='c')
        if args.freq_dstr != 'unif':
            haplotype_freqs = frequencies(args.freq_dstr, num_haplotypes,
                                          args.ratio, args.dirichlet_freqs)
            aux = np.repeat(haplotype_seqs_array, np.round(
                haplotype_freqs * 100).astype(int), axis=0)
            consensus = mode(aux, nan_policy='omit')
        else:
            consensus = mode(haplotype_seqs_array, nan_policy='omit')
        if np.any(consensus[1] < 1):
            print("At some loci the consensus base is ambiguous")
        haplotype_master_array = consensus[0][0]

    haplotype_freqs = frequencies(args.freq_dstr, num_haplotypes,
                                  args.ratio, args.dirichlet_freqs)

    # missed = np.zeros(num_haplotypes)

    df_snvs_expected = true_snvs(
        haplotype_master_array, haplotype_master, haplotype_seqs_array,
        num_haplotypes, haplotype_freqs, args.long_deletions, alphabet)

    if args.only_deletions:
        # Only look at deletions: drop other entries in expected SNVs dataframe
        if args.long_deletions:
            is_deletion = df_snvs_expected["REF"].str.len() > 1
        else:
            is_deletion = df_snvs_expected["ALT"].str.startswith('-')
        df_snvs_expected = df_snvs_expected[is_deletion]

    # Keep track of SNVs that fall within targeted regions
    df_snvs["IS_CONTAINED"] = False
    df_snvs_expected["IS_CONTAINED"] = False
    if args.long_deletions:
        deletion_length = df_snvs["REF"].str.len() - 1
        is_deletion = deletion_length > 0
        # Using 0-based indexing
        start_locus = df_snvs["POS"] - 1
        start_locus[is_deletion] += 1
        end_locus = start_locus + deletion_length - 1
        # Similarly for expected SNVs (Already uses 0-based indexing)
        deletion_length_exp = df_snvs_expected["REF"].str.len() - 1
        is_deletion_exp = deletion_length_exp > 0
        start_locus_exp = df_snvs_expected["POS"].copy()
        start_locus_exp[is_deletion_exp] += 1
        end_locus_exp = start_locus_exp + deletion_length_exp - 1
    else:
        # Handle SNVs and single-nucleotide deletions
        # Using 0-based indexing
        start_locus = df_snvs.index.get_level_values("POS") - 1
        end_locus = None
        # Similarly for expected SNVs (Already uses 0-based indexing)
        start_locus_exp = df_snvs_expected["POS"]
        end_locus_exp = None

    if args.coverage_intervals is not None:
        with open(args.coverage_intervals, 'r') as infile:
            for line in infile:
                record = line.rstrip().split('\t')
                if record[0] == args.sampleID:
                    if len(record) == 1:
                        print("Empty target region")
                        with open(args.outfile, 'w') as outfile:
                            outfile.write('ID\tTP\tFP\tFN\tTN\n')
                        return
                    regions = record[1]
                    break
        regions = regions.split(',')
        idxs = np.zeros(reference_len, dtype=bool)
        print("Reporting using 1-based indexing (and closed intervals)")
        num_loci = 0
        for r in regions:
            aux = r.split(':')
            ref_name = aux[0]
            if args.haplotype_master is not None:
                assert header == ref_name, (
                    f"Name of the reference, {ref_name}, does not agree with "
                    f"fasta file, {header}")
            aux = aux[1].split('-')
            start = int(aux[0])
            end = int(aux[1])
            if args.snv_caller == 'lofreq' or args.no_expansion:
                # Region is interpreted as a closed interval and using 1-based
                # indexing
                start -= 1
                start = max(0, start)
            elif args.snv_caller == 'shorah':
                # ShoRAH was used for SNV calling
                # Assuming 3 windows were used for SNV calling, identify
                # region that is covered by at least 2 windows (below, using
                # 0-based indexing and closed intervals)
                start_ = max(0, start - args.window_len - 1)
                end_ = min(reference_len, end + args.window_len)
                num_windows = np.floor(
                    (end_ - (start_ + args.window_len - 1)) /
                    (args.window_len // args.window_shift)) + 1
                offset = ((args.window_shift - 1) * args.window_len /
                          args.window_shift)

                start = max(0, start - offset - 1)
                # In order to identify the region which is covered by at least
                # two windows, add to the end of the first window the
                # increment multiply by the number of windows - 2 (i.e.,
                # discarding last window). In this case assuming half-open
                # interval [start, end)
                end = min(
                    reference_len, start_ + args.window_len +
                    (num_windows - 2) * (args.window_len // args.window_shift))
            # idxs[range(int(start), int(end))] = True
            # loci_region = loci[int(start):int(end)]

            # if DBG:
            #     print(f"DBG loci_true[i]: {loci_true[i]}")
            #     print(f"DBG loci_region[0]: {loci_region[0]}")
            # Here, loci are reported using 1-based indexing and a closed
            # interval
            num_loci += (end - start)
            start = int(start)
            end = int(end)
            print(f"Region with enough support: {start + 1}-{end}")

            # Mark reported and expected SNVs within the region
            is_contained = target_snvs(start, end, start_locus,
                                       args.long_deletions, end_locus)
            df_snvs["IS_CONTAINED"] = (df_snvs["IS_CONTAINED"] | is_contained)
            is_contained = target_snvs(start, end, start_locus_exp,
                                       args.long_deletions, end_locus_exp)
            df_snvs_expected["IS_CONTAINED"] = (
                df_snvs_expected["IS_CONTAINED"] | is_contained)

    else:
        loci = np.arange(reference_len)
        if args.snv_caller == 'shorah':
            idxs = np.zeros(reference_len, dtype=bool)
            offset = (args.window_len // args.window_shift)
            # Parse coverage intervals from ShoRAH output
            with open(args.coverage, 'r') as infile:
                # Look for regions at least covered by two windows
                start_w = 1
                end_w = 1
                for count, line in enumerate(infile):
                    record = line.rstrip().split("\t")
                    if count == 0:
                        start_w = int(record[2])
                        end_w = int(record[3])
                    else:
                        if int(record[2]) == start_w + offset:
                            start_w = int(record[2])
                            idxs[(start_w - 1):end_w] = True
                        else:
                            start_w = int(record[2])
                        end_w = int(record[3])

            loci_region = np.extract(idxs, loci)

        else:
            if args.coverage is not None:
                with open(args.coverage, 'r') as infile:
                    header = infile.readline().rstrip().split("\t")
                sampleID_idx = [
                    idx for idx, name in enumerate(header)
                    if args.sampleID in name
                ]
                coverage = np.loadtxt(args.coverage, dtype=int, delimiter='\t',
                                      skiprows=1, usecols=(sampleID_idx[0],))
                assert coverage.size == reference_len, (
                    "Coverage file and reference file do not have the same "
                    "number of loci")
                # Do not account for position with zero coverage for reporting
                # TP, FP, FN, and specially TN
                mask = coverage <= 0
                loci_region = loci[~mask]
            else:
                raise IOError(
                    "Expected coverage file as input when target region is not specified"
                )

        num_loci = loci_region.size
        regions = consecutive(loci_region)
        start = [el[0] for el in regions]
        end = [el[-1] + 1 for el in regions]
        for si, ei in zip(start, end):
            # Mark reported and expected SNVs within the region
            is_contained = target_snvs(si, ei, start_locus,
                                       args.long_deletions, end_locus)
            df_snvs["IS_CONTAINED"] = (df_snvs["IS_CONTAINED"] | is_contained)
            is_contained = target_snvs(si, ei, start_locus_exp,
                                       args.long_deletions, end_locus_exp)
            df_snvs_expected["IS_CONTAINED"] = (
                df_snvs_expected["IS_CONTAINED"] | is_contained)

    # Drop SNVs that fall outside of the targeted regions. Otherwise, these
    # rows will be counted toward false positives/negatives.
    df_snvs = df_snvs[df_snvs["IS_CONTAINED"]]
    df_snvs_expected = df_snvs_expected[df_snvs_expected["IS_CONTAINED"]]

    if args.output_true:
        output_file = os.path.join(outdir, 'true_snvs.tsv')
        # Report using 1-based indexing
        df_snvs_expected["POS"] += 1
        df_snvs_expected.to_csv(
            output_file, sep="\t",
            columns=["POS", "REF", "ALT", "FREQ", "HAPLOTYPES"],
            header=["Loci", "Reference", "Variant", "Frequency", "Haplotypes"],
            index=False, compression=None)

    # join on POS and ALT
    df_pairs = df_snvs_expected.merge(
        df_snvs, how="outer", on=["POS", "ALT", "REF"],
        suffixes=["_exp", "_rep"])

    FN_mask = df_pairs["INFO"].isnull()
    FN = sum(FN_mask)

    FP_mask = df_pairs["HAPLOTYPES"].isnull()
    FP = sum(FP_mask)

    TP_mask = ~FN_mask & ~FP_mask
    TP = sum(TP_mask)

    TN = num_loci - len(df_pairs["POS"].value_counts())
    # Sensitivity
    if TP or FN:
        print("Sensitivity: {:.6f}".format(TP / (TP + FN)))

    # Precision
    if TP or FP:
        print("Precision: {:.6f}".format(TP / (TP + FP)))

    # Specificity
    if TN or FP:
        print("Specificity: {:.6f}".format(TN / (TN + FP)))

    print("TP: ", TP)
    print("FP: ", FP)
    print("FN: ", FN)
    print("TN: ", int(TN))
    # print("Number of FN per haplotype: ", missed)

    # Write to output file
    with open(args.outfile, 'w') as outfile:
        outfile.write('ID\tTP\tFP\tFN\tTN\n')
        outfile.write(f'{args.sampleID}\t{TP}\t{FP}\t{FN}\t{int(TN)}\n')

    # output_file = os.path.join(outdir, 'FN_per_haplotype.tsv')
    # with open(output_file, 'w') as outfile:
    #     for idx, name in enumerate(haplotype_ids):
    #         aux = name.split(' ')[0]
    #         outfile.write(f'{aux}\t{missed[idx]}\n')

    output_file = os.path.join(outdir, 'TP_frequencies.tsv')
    df_pairs[TP_mask].to_csv(output_file, sep="\t",
                             columns=["POS", "REF", "ALT", "FREQ_exp",
                                      "FREQ_rep", "INFO"],
                             header=["Loci", "Reference", "Variant",
                                     "Freq (expected)", "Freq (reported)",
                                     "Info"],
                             index=False, compression=None)

    output_file = os.path.join(outdir, 'FP_frequencies.tsv')
    df_pairs[FP_mask].to_csv(output_file, sep="\t",
                             columns=["POS", "REF", "ALT", "FREQ_rep", "INFO"],
                             header=["Loci", "Reference", "Variant", "Freq",
                                     "Info"],
                             index=False, compression=None)

    output_file = os.path.join(outdir, 'FN_frequencies.tsv')
    df_pairs[FN_mask].to_csv(output_file, sep="\t",
                             columns=["POS", "REF", "ALT", "FREQ_exp",
                                      "HAPLOTYPES"],
                             header=["Loci", "Reference", "Variant", "Freq",
                                     "Haplotypes"],
                             index=False, compression=None)


if __name__ == '__main__':
    main()
