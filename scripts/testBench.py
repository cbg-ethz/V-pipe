#!/usr/bin/env python3

import os
import argparse
import csv
from math import trunc
from coverageIntervals import read_fasta

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sh
import numpy as np


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
        "--no-shorah", required=False, default=False, action='store_true',
        dest='snv_caller',
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


def true_snvs(haplotype_master, haplotype_seqs, num_haplotypes, haplotype_freqs,
              alphabet):
    """
    Extract expected SNVs using the MSA of the true haplotype sequences and
    the reference sequence
    """
    loci = np.arange(haplotype_master.size)
    haplotype_idx = np.arange(num_haplotypes)
    variants = haplotype_master != haplotype_seqs

    loci_true = []
    ref_true = []
    snvs_true = []
    freq_true = []
    haps_true = []
    for idx in loci:
        idxs = variants[:, idx]
        if np.any(idxs):
            var = haplotype_seqs[idxs, idx]
            snv_freq = haplotype_freqs[idxs]
            haplotypeIDs = haplotype_idx[idxs]
            if np.sum(idxs) == 1:
                loci_true.append(loci[idx])
                ref_true.append(haplotype_master[idx].decode())
                snvs_true.append(var[0].decode())
                freq_true.append(snv_freq[0])
                haps_true.append(str(haplotypeIDs[0]))
            else:
                for a in alphabet:
                    idxs_a = var == a
                    hapIDs = haplotypeIDs[idxs_a]
                    if np.sum(idxs_a) > 0:
                        loci_true.append(loci[idx])
                        ref_true.append(haplotype_master[idx].decode())
                        snvs_true.append(a)
                        freq_true.append(np.sum(snv_freq[idxs_a]))
                        haps_true.append(
                            str(hapIDs[0]) if len(idxs_a) == 1 else ','.join((
                                str(x) for x in hapIDs)))

    snvs_true = np.array(snvs_true, dtype='c')

    return loci_true, ref_true, snvs_true, freq_true, haps_true


def inferred_snvs(snvs_file):
    "Compute average frequency for SNVs called using ShoRAH"
    loci_inferred = []
    ref_inferred = []
    snvs_inferred = []
    freq_inferred = []
    loci_inferred_tmp = []
    ref_inferred_tmp = []
    snvs_inferred_tmp = []
    freq_inferred_tmp = []
    del_len = 0

    extension = os.path.splitext(snvs_file)[1][1:]
    if extension == "vcf":
        with open(snvs_file) as infile:
            for line in infile:
                record = line.rstrip()
                if record and record[0] != '#':
                    row = record.split("\t")
                    freq = row[7]
                    freq = freq.split("AF=")[1]
                    freq = freq.split(";")[0]
                    # first base corresponds to the locus after which the
                    # deletion starts
                    ref_len = len(row[3]) - 1
                    var_len = len(row[4])

                    # Check that variant do not correspond to a insertion
                    if var_len == 1:
                        # Check if variant corresponds to a deletion
                        if ref_len >= 1:
                            # first base corresponds to the locus after which
                            # the deletion starts
                            locus = int(row[1])
                            # report deletions per loci w.r.t reference
                            if del_len > 0:
                                idx = locus - loci_inferred_tmp[0]
                                assert idx >= 0
                                loci_inferred = loci_inferred + \
                                    loci_inferred_tmp[:idx]
                                ref_inferred = ref_inferred + \
                                    ref_inferred_tmp[:idx]
                                snvs_inferred = snvs_inferred + \
                                    snvs_inferred_tmp[:idx]
                                freq_inferred = freq_inferred + \
                                    freq_inferred_tmp[:idx]
                                # shorten the current deletion
                                loci_inferred_tmp = loci_inferred_tmp[idx:]
                                ref_inferred_tmp = ref_inferred_tmp[idx:]
                                snvs_inferred_tmp = snvs_inferred_tmp[idx:]
                                freq_inferred_tmp = freq_inferred_tmp[idx:]
                                del_len = len(loci_inferred_tmp)
                            # There are two options:
                            # (1) the second deletion outruns the first one, or
                            # (2) the second deletion is fully contained in the
                            #     first one
                            if ref_len > del_len:
                                # Add deleted bases to current deletion
                                for idx in range(ref_len):
                                    if idx >= del_len:
                                        locus = int(row[1]) + idx
                                        loci_inferred_tmp.append(locus)
                                        # ALL deleted bases w.r.t. reference
                                        # are reported. However, the first base
                                        # corresponds to the locus after which
                                        # the deletion starts
                                        ref_inferred_tmp.append(row[3][idx + 1])
                                        snvs_inferred_tmp.append('-')
                                        freq_inferred_tmp.append(float(freq))
                                    else:
                                        freq_inferred_tmp[idx] += float(freq)
                            else:
                                # Add frequencies
                                for idx in range(ref_len):
                                    freq_inferred_tmp[idx] += float(freq)
                            del_len = len(loci_inferred_tmp)
                        else:
                            # Lofreq uses 1-based indexing
                            locus = int(row[1]) - 1
                            # check if a deletion is to be inserted before this
                            # locus
                            if del_len > 0:
                                # There are two options:
                                # (1) The next locus happens after the last
                                #     locus of the deletion, or
                                # (2) in between the deletion
                                if loci_inferred_tmp[-1] <= locus:
                                    loci_inferred = loci_inferred + loci_inferred_tmp
                                    ref_inferred = ref_inferred + ref_inferred_tmp
                                    snvs_inferred = snvs_inferred + snvs_inferred_tmp
                                    freq_inferred = freq_inferred + freq_inferred_tmp
                                    loci_inferred_tmp = []
                                    ref_inferred_tmp = []
                                    snvs_inferred_tmp = []
                                    freq_inferred_tmp = []
                                else:
                                    # insert leading loci - loci before snv
                                    idx = locus - loci_inferred_tmp[0] + 1
                                    # Deletions are reported at the preceding
                                    # locus, so it can be that an snv at this
                                    # locus is reported after the deletion
                                    if idx > 0:
                                        loci_inferred = loci_inferred + \
                                            loci_inferred_tmp[:idx]
                                        ref_inferred = ref_inferred + \
                                            ref_inferred_tmp[:idx]
                                        snvs_inferred = snvs_inferred + \
                                            snvs_inferred_tmp[:idx]
                                        freq_inferred = freq_inferred + \
                                            freq_inferred_tmp[:idx]
                                        loci_inferred_tmp = loci_inferred_tmp[idx:]
                                        ref_inferred_tmp = ref_inferred_tmp[idx:]
                                        snvs_inferred_tmp = snvs_inferred_tmp[idx:]
                                        freq_inferred_tmp = freq_inferred_tmp[idx:]
                            # insert the snv
                            loci_inferred.append(locus)
                            ref_inferred.append(row[3])
                            snvs_inferred.append(row[4])
                            freq_inferred.append(float(freq))
                            del_len = len(loci_inferred_tmp)
    else:
        with open(snvs_file, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            next(reader, None)    # skip the header
            for row in reader:
                locus = int(row[1]) - 1    # ShoRAH uses 1-based indexing
                freq = np.zeros(shape=3)
                for idx in np.arange(3):
                    aux = row[idx + 4]
                    if aux == '-':
                        freq[idx] = 0.0
                    elif aux == '*':
                        freq[idx] = np.nan
                    else:
                        freq[idx] = float(aux)
                freq = np.nanmean(freq)
                loci_inferred.append(locus)
                ref_inferred.append(row[2])
                snvs_inferred.append(row[3])
                freq_inferred.append(freq)

    snvs_inferred = np.array(snvs_inferred, dtype='c')

    return loci_inferred, ref_inferred, snvs_inferred, freq_inferred


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


def get_FN(locus, i, i_max, loci_true, snvs_true, freq_true, FN, FN_freq,
           missed, haps_true):
    "Count false negatives"
    while loci_true[i] == locus:
        FN += 1
        print("At locus {}, SNV missed {} ({:.6f})".format(
            loci_true[i] + 1, snvs_true[i].decode('UTF-8'), freq_true[i]))
        FN_freq = FN_freq + \
            [[loci_true[i] + 1,
                snvs_true[i].decode('UTF-8'), freq_true[i], haps_true[i]]]
        missed[[int(x) for x in haps_true[i].split(',')]] += 1
        i += 1
        if i == i_max:
            break
    return FN, FN_freq, i, missed


def get_FN_TN(locus, i, i_max, loci_true, snvs_true, freq_true, FN, TN, FN_freq,
              missed, haps_true):
    """
    Handle cases in which there are no more reported SNVs, but some are
    expected
    """
    if loci_true[i] == locus:
        FN, FN_freq, i, missed = get_FN(
            locus, i, i_max, loci_true, snvs_true, freq_true, FN, FN_freq,
            missed, haps_true)
    else:
        if DBG:
            print(f"At locus {locus + 1}, true negative reported")
        TN += 1
    return FN, TN, FN_freq, i, missed


def get_FP(locus, j, j_max, loci_inferred, snvs_inferred, freq_inferred, FP,
           FP_freq):
    "Count false positives"
    while loci_inferred[j] == locus:
        FP += 1
        print(
            "At locus {}, technical error / artefact reported as true variant: {} ({:.6f})"
            .format(loci_inferred[j] + 1, snvs_inferred[j].decode('UTF-8'),
                    freq_inferred[j]))
        FP_freq = FP_freq + [[loci_inferred[j] + 1,
                              snvs_inferred[j].decode('UTF-8'),
                              freq_inferred[j]]]
        j += 1
        if j == j_max:
            break
    return FP, FP_freq, j


def get_FP_TN(locus, j, j_max, loci_inferred, snvs_inferred, freq_inferred, FP,
              TN, FP_freq):
    """
    Handle cases in which there are no more expected SNVs, but nevertheless
    reported by the caller
    """
    if loci_inferred[j] == locus:
        FP, FP_freq, j = get_FP(locus, j, j_max, loci_inferred, snvs_inferred,
                                freq_inferred, FP, FP_freq)
    else:
        if DBG:
            print(f"At locus {locus + 1}, true negative reported")
        TN += 1
    return FP, TN, FP_freq, j


def consecutive(array, stepsize=1):
    return np.split(array, np.where(np.diff(array) != stepsize)[0] + 1)


def get_performance(loci_true, loci_inferred, snvs_true, snvs_inferred,
                    freq_true, freq_inferred, haps_true, num_haplotypes,
                    loci_region, i, j, TP, FP, TN, FN, TP_freq, FP_freq,
                    FN_freq, missed, coverage_file=False, regions=None):
    i_max = len(loci_true)
    j_max = len(loci_inferred)

    while loci_true[i] < loci_region[0]:
        i += 1
    if DBG:
        print(f"DBG loci_true[i]: {loci_true[i]}")
    # Inferred loci can be outside of the target region, e.g., when
    # reporting metrics based on intersection between tools
    while j < j_max and loci_inferred[j] < loci_region[0]:
        j += 1

    idx_region = 0
    for idx in loci_region:

        if coverage_file:
            if len(regions) > idx_region:
                if idx == regions[idx_region][0]:
                    while loci_true[i] < idx:
                        i += 1
                    print("Region with enough support: {:d}-{:d}".format(
                        int(regions[idx_region][0]) + 1,
                        int(regions[idx_region][-1] + 1)))
                    idx_region += 1

        if i == i_max or j == j_max:
            if i == i_max and j < j_max:
                # There are no more expected SNVs, but nevertheless reported
                # by the caller
                FP, TN, FP_freq, j = get_FP_TN(
                    idx, j, j_max, loci_inferred, snvs_inferred, freq_inferred,
                    FP, TN, FP_freq)

            if j == j_max and i < i_max:
                # There are no more reported SNVs, but some are expected
                FN, TN, FN_freq, i, missed = get_FN_TN(
                    idx, i, i_max, loci_true, snvs_true, freq_true, FN, TN,
                    FN_freq, missed, haps_true)

            if i == i_max and j == j_max:
                break
        else:
            assert loci_true[i] >= idx
            assert loci_inferred[j] >= idx

            if loci_true[i] == idx and loci_inferred[j] == idx:
                while loci_true[i] == idx and loci_inferred[j] == idx:
                    if snvs_true[i] == snvs_inferred[j]:
                        TP += 1
                        if DBG:
                            print(
                                "At locus {}, true positive reported: {} ({:.6f})"
                                .format(loci_inferred[j] + 1,
                                        snvs_inferred[j].decode('UTF-8'),
                                        freq_inferred[j]))
                        TP_freq = TP_freq + \
                            [[loci_inferred[j] + 1,
                                snvs_inferred[j].decode('UTF-8'),
                                freq_true[i], freq_inferred[j]]]
                        i += 1
                        j += 1
                    elif snvs_true[i] > snvs_inferred[j]:
                        FP += 1
                        print(
                            "At locus {}, technical error / artefact reported as true variant: {} ({:.6f})"
                            .format(loci_inferred[j] + 1,
                                    snvs_inferred[j].decode('UTF-8'),
                                    freq_inferred[j]))
                        FP_freq = FP_freq + \
                            [[loci_inferred[j] + 1,
                                snvs_inferred[j].decode('UTF-8'),
                                freq_inferred[j]]]
                        j += 1
                    elif snvs_true[i] < snvs_inferred[j]:
                        FN += 1
                        print("At locus {}, SNV missed {} ({:.6f})".format(
                            loci_true[i] + 1, snvs_true[i].decode('UTF-8'),
                            freq_true[i]))
                        missed[[int(x) for x in haps_true[i].split(',')]] += 1
                        FN_freq = FN_freq + \
                            [[loci_true[i] + 1,
                                snvs_true[i].decode('UTF-8'), freq_true[i],
                                haps_true[i]]]
                        i += 1
                    if i == i_max or j == j_max:
                        break

                if i < i_max:
                    FN, FN_freq, i, missed = get_FN(
                        idx, i, i_max, loci_true, snvs_true, freq_true, FN,
                        FN_freq, missed, haps_true)
                if j < j_max:
                    FP, FP_freq, j = get_FP(
                        idx, j, j_max, loci_inferred, snvs_inferred,
                        freq_inferred, FP, FP_freq)

                if i == i_max and j == j_max:
                    # There are no expected nor reported SNVs. Hence, the rest
                    # of the positions correspond to true negatives
                    break
            elif loci_true[i] == idx:
                FN, FN_freq, i, missed = get_FN(
                    idx, i, i_max, loci_true, snvs_true, freq_true, FN, FN_freq,
                    missed, haps_true)
            elif loci_inferred[j] == idx:
                FP, FP_freq, j = get_FP(
                    idx, j, j_max, loci_inferred, snvs_inferred, freq_inferred,
                    FP, FP_freq)
            else:
                if DBG:
                    print(f"At locus {idx + 1}, true negative reported")
                TN += 1

    # Add positions that were not reported as polymorphic by the caller and
    # were not expected as true SNVs
    if DBG:
        if loci_region[-1] > idx:
            print(f"Loci {idx + 1}-{loci_region[-1]}, reported as true negative")
    TN += loci_region[-1] - idx

    return TP, FP, TN, FN, TP_freq, FP_freq, FN_freq, missed, i, j


def main():

    args = parse_args()

    alphabet = ['-', 'A', 'C', 'G', 'T']
    alphabet = np.array(alphabet, dtype='c')

    # Compute average frequency for SNVs called using ShoRAH
    loci_inferred, ref_inferred, snvs_inferred, freq_inferred = inferred_snvs(
        args.snvs)
    if not loci_inferred:
        print("No called SNVs")
        with open(args.outfile, 'w') as outfile:
            outfile.write('ID\tTP\tFP\tFN\tTN\n')
        return
    
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
                cohort_consensus = SeqRecord(Seq(''.join(haplotype_master_array)),
                                             id=header, description="")
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
    # True haplotypes - expected SNVs
    loci_true, ref_true, snvs_true, freq_true, haps_true = true_snvs(
        haplotype_master_array, haplotype_seqs_array, num_haplotypes,
        haplotype_freqs, alphabet)

    if args.output_true:
        output_file = os.path.join(outdir, 'true_snvs.tsv')
        with open(output_file, 'w') as outfile:
            outfile.write('Locus\tRef\tVar\tFreq\tHaplotypes\n')
            for idx in range(len(loci_true)):
                outfile.write('{}\t{}\t{}\t{}\t{}\n'.format(
                    loci_true[idx] + 1, ref_true[idx],
                    snvs_true[idx].decode('utf-8'), freq_true[idx],
                    haps_true[idx]))

    missed = np.zeros(num_haplotypes)
    # TP: loci that are truly polymorphic
    TP = 0
    # FP: technical error reported as SNVs
    FP = 0
    # TN: loci that are not polymorphic
    TN = 0
    # FN: SNVs that are missed
    FN = 0
    # SNV frequencies
    TP_freq = []
    FP_freq = []
    FN_freq = []

    loci = np.arange(reference_len)
    i = 0
    j = 0

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
            if args.snv_caller:
                # Region is interpreted as a closed interval and using 1-based
                # indexing
                start -= 1
                start = max(0, start)
            else:
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
            idxs[range(int(start), int(end))] = True
            loci_region = loci[int(start):int(end)]

            if DBG:
                print(f"DBG loci_true[i]: {loci_true[i]}")
                print(f"DBG loci_region[0]: {loci_region[0]}")
            # Here, loci are reported using 1-based indexing and a closed
            # interval
            print("Region with enough support: {:d}-{:d}".format(
                trunc(start) + 1, trunc(end)))

            TP, FP, TN, FN, TP_freq, FP_freq, FN_freq, missed, i, j = get_performance(
                loci_true, loci_inferred, snvs_true, snvs_inferred, freq_true,
                freq_inferred, haps_true, num_haplotypes, loci_region, i, j,
                TP, FP, TN, FN, TP_freq, FP_freq, FN_freq, missed)

        loci = loci[idxs]
        if loci_inferred[0] < loci[0] or loci_inferred[-1] > loci[-1]:
            print("Warning: some reported SNVs are outside the target region."
                  " It can happen when target region is smaller than region"
                  " where SNVs were called.")
    else:
        if not args.snv_caller:
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

        regions = consecutive(loci_region)

        TP, FP, TN, FN, TP_freq, FP_freq, FN_freq, missed, i, j = get_performance(
            loci_true, loci_inferred, snvs_true, snvs_inferred, freq_true,
            freq_inferred, haps_true, num_haplotypes, loci_region, i, j, TP,
            FP, TN, FN, TP_freq, FP_freq, FN_freq, missed,
            coverage_file=True, regions=regions)

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
    print("TN: ", TN)
    print("Number of FN per haplotype: ", missed)

    # Write to output file
    with open(args.outfile, 'w') as outfile:
        outfile.write('ID\tTP\tFP\tFN\tTN\n')
        outfile.write(f'{args.sampleID}\t{TP}\t{FP}\t{FN}\t{TN}\n')

    output_file = os.path.join(outdir, 'FN_per_haplotype.tsv')
    with open(output_file, 'w') as outfile:
        for idx, name in enumerate(haplotype_ids):
            aux = name.split(' ')[0]
            outfile.write(f'{aux}\t{missed[idx]}\n')

    output_file = os.path.join(outdir, 'TP_frequencies.tsv')
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerow(['Loci', 'Variant', 'Freq', 'Inferred freq'])
        writer.writerows(TP_freq)

    output_file = os.path.join(outdir, 'FP_frequencies.tsv')
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerow(['Loci', 'Variant', 'Inferred freq'])
        writer.writerows(FP_freq)

    output_file = os.path.join(outdir, 'FN_frequencies.tsv')
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerow(['Loci', 'Variant', 'Freq'])
        writer.writerows(FN_freq)


if __name__ == '__main__':
    main()
