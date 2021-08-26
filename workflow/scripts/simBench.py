#!/usr/bin/env python3

import os
import sys
import argparse
import collections
import shutil
from itertools import zip_longest
from alignmentIntervals import read_fasta

import sh
import numpy as np

__author__ = "Susana Posada-Cespedes"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"


class Tree(object):
    def __init__(self, sequence, max_depth, parent=None):
        self.sequence = sequence
        self.max_depth = max_depth
        self.parent = parent
        self.children = []

    @property
    def depth(self):
        return 0 if not self.parent else 1 + self.parent.depth

    def get_leaves(self):
        leaves = []

        def _get_leaf(node):
            if node:
                if len(node.children) == 0:
                    leaves.append(node.sequence)
                for child in node.children:
                    _get_leaf(child)

        _get_leaf(self)

        return leaves


def parse_args():
    """ Set up the parsing of command-line arguments """

    parser = argparse.ArgumentParser(
            description="Script for simulating data",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "-f", required=False, default=None, metavar='FASTA',
        dest='haplotype_seqs',
        help="Fasta file containing sequences of the master haplotype(s)"
    )
    parser.add_argument(
        "-u", required=False, default=False, action='store_true',
        dest='use_master',
        help="Generate haplotype sequences from master sequence"
    )
    parser.add_argument(
        "-g", required=False, default=1000, metavar='INT',
        dest='genome_length', type=int,
        help="Length of the genome or target region"
    )
    parser.add_argument(
        "-n", required=False, default=5, metavar='INT',
        dest="num_haplotypes", type=int, help="Number of underlying haplotypes"
    )
    parser.add_argument(
        "-c", required=False, default='500', metavar='INT,INT,...',
        dest='coverage', type=str,
        help="Average coverage per locus. When more than one value is given, "
             "it is interpreted as the coverage per haplotype"
    )
    parser.add_argument(
        "-nr", required=False, default=False, action='store_true',
        dest='num_reads',
        help="Indicate if argument -c should be interpreted as number of reads"
    )
    parser.add_argument(
        "-l", required=False, default=250, metavar='INT', dest='read_length',
        type=int, help="Length (in base pairs) of simulated reads"
    )
    parser.add_argument(
        "-p", required=False, default=False, action='store_true',
        dest='paired', help="Indicate whether to simulate paired-end reads"
    )
    parser.add_argument(
        "-m", required=False, default='600,100', metavar='mean,sd',
        dest='fragment_size',
        help="Average fragment size and its standard deviation"
    )
    parser.add_argument(
        "-d", required=False, default='unif', metavar='STR', dest='freq_dstr',
        type=str, choices=['unif', 'geom', 'dirichlet', 'cust'],
        help="Distribution of haplotype frequencies"
    )
    parser.add_argument(
        "-gr", required=False, default=0.75, metavar='FLOAT',
        dest='geom_ratio', type=float,
        help="Sucess probability for the geometric distribution"
    )
    parser.add_argument(
        "-dc", required=False, default=None, metavar='FLOAT,FLOAT,...',
        dest='dirichlet_alpha', type=str,
        help="Concentration parameters of the Dirichlet distribution. For "
             "values > 1, more evenly distributed distributions are "
             "generated. For values < 1, the vast majority of the mass will "
             "be concentrated in a few of the values."
    )
    parser.add_argument(
        "-mr", required=False, default=0.1, metavar='FLOAT',
        dest='mutation_rate', type=float,
        help="Mutation rate for generating variants/haplotypes. If greater "
             "than 1 is interpreted as the desired number of mutations to "
             "introduce per sequence"
    )
    parser.add_argument(
        "-dr", required=False, default=0.0, metavar='FLOAT',
        dest='deletion_rate', type=float,
        help="Deletion rate for generating variants/haplotypes. If greater "
             "than 1 is interpreted as the desired number of deletions to "
             "introduce per sequence"
    )
    parser.add_argument(
        "-ir", required=False, default=0.0, metavar='FLOAT',
        dest='insertion_rate', type=float,
        help="Insertion rate for generating variants/haplotypes. If greater "
             "than 1 is interpreted as the desired number of insertions to "
             "introduce per sequence")
    parser.add_argument(
        "-fr", required=False, default=False, action='store_true',
        dest='no_FR',
        help="Indicate if frame-shift mutations should be avoided, and rather "
             "long deletions introduced"
    )
    parser.add_argument(
        "-dl", required=False, default=None, metavar='INT',
        dest='deletion_length', type=int, help="Deletion length in bp"
    )
    parser.add_argument(
        "-q", required=False, default=False, action='store_true',
        dest='quality',
        help="Indicate whether to simulate high-quality reads"
    )
    parser.add_argument(
        "-art", required=False, default="art_illumina", metavar='PATH',
        dest='art', type=str,
        help="Path to binaries for the read simulator ART"
    )
    parser.add_argument(
        "--tree-like", required=False, default=False, action='store_true',
        dest='tree_like', help="Indicate whether to generate haplotype by "
                               "emulating a tree-like relationship among the "
                               "sequences"
    )
    parser.add_argument(
        "-s", required=False, default=None, metavar='INT', dest='seed',
        type=int, help="Set seed for reproducibility"
    )
    parser.add_argument(
        "-v", required=False, default=False, action='store_true',
        dest='verbose', help="Verbose output"
    )
    parser.add_argument(
        "-oh", required=False, default=None, metavar='DIR',
        dest='outdir_haps', type=str,
        help="Output directory for haplotype sequences"
    )
    parser.add_argument(
        "-or", required=False, default=None, metavar='DIR',
        dest='outdir_reads', type=str, help="Output directory for reads"
    )
    parser.add_argument(
        "-o", required=False, default="all", metavar='STR', dest='output',
        type=str, choices=['master', 'haplotypes', 'reads', 'all'],
        help="Specify the output: 'master' for generating the master sequence, "
             "'haplotypes' for simulating haplotypes, 'reads' for simulating "
             "reads, and 'all' for three outputs. Only specify 'reads' if input"
             " files were generated before with the same script"
    )

    return parser.parse_args()


def grouper(n, iterable, fillvalue=None):
    """
    Collect data into fixed-length chunks or blocks
    Source: itertools recipes
    """
    # grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def mutations(num_mutations, seed=None):
    """
    Generate point mutations as integers, and translate integers into
    nucleotides using the mapping defined by:
    index base
    0     A
    1     C
    2     G
    3     T
    """
    # Set seed for reproducibility
    np.random.seed(seed)

    int_seq = np.random.randint(0, 3, size=num_mutations)

    alphabet_idxs = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    nt_seq = np.zeros(int_seq.size, dtype='c')
    for idx, nt in alphabet_idxs.items():
        nt_seq[int_seq == idx] = nt
    return nt_seq


# Adapted from: http://hplgit.github.io/bioinf-py/doc/pub/bioinf-py.html
def mutate(haplotype, mutation_rate, del_rate, ins_rate, noFR=True,
           del_len=None, seed=None, verbose=False):
    """
    Mutate haplotype into a related sequence by introducing mutations at rate
    mutation_rate and deletions at rate del_rate.
    """
    seed = seed if seed is not None else np.random.randint(1, 1000, size=1)[0]
    # Set seed for reproducibility
    np.random.seed(seed)

    haplotype = np.array(haplotype, dtype='c')
    if mutation_rate >= 1:
        num_mutations = int(mutation_rate)
    else:
        num_mutations = int(mutation_rate * haplotype.size)
    if del_rate >= 1:
        num_deletions = int(del_rate)
    else:
        num_deletions = int(del_rate * haplotype.size)
    if ins_rate >= 1:
        num_insertions = int(ins_rate)
    else:
        num_insertions = int(ins_rate * haplotype.size)

    if verbose:
        print("Generating mutated sequence ...")

    # Draw positions at random where to introduce point mutations
    mutated_loci = np.random.randint(0, haplotype.size - 1, size=num_mutations)
    # Draw the new bases.
    mutated_nt = mutations(num_mutations, seed=seed + 11)
    counter = 1
    while np.any(haplotype[mutated_loci] == mutated_nt):
        if verbose:
            print("Some proposed mutations match the reference")
        aux_nt = mutations(
            np.sum(haplotype[mutated_loci] == mutated_nt), seed=seed + counter)
        mutated_nt[haplotype[mutated_loci] == mutated_nt] = aux_nt
        counter += 1

    if verbose:
        ref_nt = haplotype[mutated_loci]
        for idx in range(num_mutations):
            print(
                f"{mutated_loci[idx]}:{ref_nt[idx].decode()}->"
                f"{mutated_nt[idx].decode()}")
    haplotype[mutated_loci] = mutated_nt

    if noFR or del_len:
        # Introduce large indels (multiple of 3, such that they don't introduce
        # frame-shift mutations)
        for i in range(num_deletions):
            if del_len is None:
                del_len_cur = np.random.choice([3, 6, 9, 12], 1)[0]
            else:
                del_len_cur = del_len
            del_locus = np.random.randint(
                0, haplotype.size - del_len_cur, size=1)[0]
            if verbose:
                del_out = haplotype[del_locus:del_locus + del_len_cur].tostring().decode()
                print(f"{del_locus}:del-{del_out}")
            haplotype[del_locus:del_locus + del_len_cur] = '-'

        for i in range(num_insertions):
            ins_locus = np.random.randint(
                0, haplotype.size - 1, size=1)
            ins_len = np.random.choice([3, 6, 9, 12], 1)
            ins_seq = mutations(ins_len[0], seed=seed + 31)
            aux = np.append(haplotype[:ins_locus[0]], ins_seq)
            haplotype = np.append(aux, haplotype[ins_locus[0]:])
            if verbose:
                # print the locus (0-indexing) after which the insertion is
                # placed and the sequence
                print(
                    f"{ins_locus[0] - 1}:ins-{ins_seq.tostring().decode()}")
    else:
        # Draw positions at random where to introduce deletions
        del_loci = np.random.randint(
            0, haplotype.size - 1, size=num_deletions)
        if verbose:
            ref_nt = haplotype[del_loci]
            for idx in range(num_deletions):
                print(f"{del_loci[idx]}:del-{ref_nt[idx].decode()}")
        haplotype[del_loci] = '-'

        for i in range(num_insertions):
            ins_locus = np.random.randint(
                0, haplotype.size - 1, size=1)
            ins_seq = mutations(1, seed=seed + 41)
            aux = np.append(haplotype[:ins_locus[0]], ins_seq)
            haplotype = np.append(aux, haplotype[ins_locus[0]:])
            if verbose:
                # print the locus (0-indexing) after which the insertion is
                # placed and the sequence
                print(
                    f"{ins_locus[0] - 1}:ins-{ins_seq[0].decode()}")

    return haplotype.tostring().decode('utf-8')


def leaf(tree, mutation_rate, del_rate, ins_rate, noFR=True, del_len=None,
         seed=None, verbose=False):
    if tree.depth < tree.max_depth:
        sequences = [mutate(
            tree.sequence, mutation_rate, del_rate, ins_rate, noFR,
            del_len, seed + i + tree.depth, verbose) for i in range(2)]
        tree.children = [Tree(sequences[0], tree.max_depth, tree),
                         Tree(sequences[1], tree.max_depth, tree)]
        for child in tree.children:
            leaf(child, mutation_rate, del_rate, ins_rate, noFR, del_len, seed,
                 verbose)
    else:
        return


def sim_haplotypes(length, weights=None, s=None):
    """
    Generate random haplotypes (DNA sequences) of a given length.
    Weights can be provided to enforce a certain composition of the sequence.
    These weights should sum up to 1, e.g., (('A', 0.1), ('C', 0.2),
    ('G', 0.4), ('T', 0.3))
    """

    # Set seed for reproducibility
    np.random.seed(s)

    if weights is None:
        new_haplotype = mutations(length, seed=s)
    else:
        assert sum(w for nt, w in weights), "Weights do not sum up to 1"
        alphabet = ''.join(nt * int(w * 100) for nt, w in weights)
        alphabet = alphabet.upper()
        alphabet = np.array(alphabet, dtype='c')
        new_haplotype = np.random.choice(alphabet, length)

    new_haplotype = new_haplotype.tostring().decode('utf-8')

    return new_haplotype


def sim_master(length, outdir_haps, seed=None):
    """Simulate master sequence and save record to FASTA file"""
    master_seq = sim_haplotypes(length=length, s=seed)
    output_file = os.path.join(outdir_haps, "haplotype_master.fasta")
    with open(output_file, 'w') as outfile:
        outfile.write(f">master\n{master_seq}\n")

    return master_seq


def write_fasta(haplotype_seqs, outdir):
    fasta_record = collections.namedtuple("fasta_record", "id seq")
    output_files = []
    for idx in range(len(haplotype_seqs)):
        haplotype_id = ''.join(("haplotype", str(idx)))
        seq = fasta_record(id=haplotype_id, seq=haplotype_seqs[idx])
        output_file = os.path.join(outdir, ''.join((haplotype_id, ".fasta")))
        output_files.append(output_file)

        with open(output_file, 'w') as outfile:
            outfile.write(">{}\n{}\n".format(seq.id, seq.seq))

    sh.cat(output_files, _out=os.path.join(outdir, "haplotypes.fasta"))


def sim_reads(art, haplotype_seq, coverage, read_len, fragment_mean,
              fragment_sd, outprefix='', paired=True, highQ=True,
              num_reads=False, seed=None, verbose=False):
    # NOTE: art_illumina additional interesting options: -na (no aln file)
    cmd = sh.Command(art)
    cmd = cmd.bake('-i', haplotype_seq)
    if num_reads:
        cmd = cmd.bake('-c', coverage)
    else:
        cmd = cmd.bake('-f', coverage)
    cmd = cmd.bake('-l', read_len)
    cmd = cmd.bake('-o', outprefix)
    cmd = cmd.bake('-ef')   # generate error free reads
    cmd = cmd.bake('-sam')  # generate SAM alignment file
    # use 'M' instead of '=/X' in the CIGAR for alignment match/mismatch
    cmd = cmd.bake('-M')
    # the MSv3 profile seems to generate poorer quality reads
    # cmd = cmd.bake('-ss', 'MSv3')

    if paired:
        cmd = cmd.bake('-p')
        cmd = cmd.bake('-m', fragment_mean)
        cmd = cmd.bake('-s', fragment_sd)

    if highQ:
        cmd = cmd.bake('-qL', 28)
        # ngshmmalign range '!' (=0) and 'I' (=41).
        cmd = cmd.bake('-qU', 40)
        cmd = cmd.bake('-qs', 5)
        cmd = cmd.bake('-qs2', 5)

    if seed is not None:
        cmd = cmd.bake('-rs', seed)

    if verbose:
        print(cmd)
    cmd()


def sed(expr, target_file, verbose=False):
    cmd = sh.sed.bake('-i')
    cmd = cmd.bake('-e', expr)
    cmd = cmd.bake(target_file)

    if verbose:
        print(cmd)
    cmd()


def print_params(mutation_rate, deletion_rate, insertion_rate):

    if mutation_rate >= 1:
        print(f"Number of mutations: {int(mutation_rate)}")
    else:
        print(f"Mutation rate: {mutation_rate}")
    if deletion_rate >= 1:
        print(f"Number of deletions: {int(deletion_rate)}")
    else:
        print(f"Deletion rate: {deletion_rate}")
    if insertion_rate >= 1:
        print(f"Number of insertions per haplotype: "
              f"{int(insertion_rate)}")
    else:
        print(f"Insertion rate: {insertion_rate}")


def main():
    args = parse_args()

    outdir_haps = args.outdir_haps if args.outdir_haps is not None else os.getcwd()
    outdir_reads = args.outdir_reads if args.outdir_reads is not None else os.getcwd()
    seed = args.seed if args.seed is not None else np.random.randint(
        1, 1000, size=1)[0]
    num_haplotypes = args.num_haplotypes

    if args.output == "master":
        # Simulate master sequence
        haplotype_seq = sim_master(
            length=args.genome_length, outdir_haps=outdir_haps, seed=seed)
    elif args.output == "haplotypes" or args.output == "all":
        # Simulate true underlying haplotypes
        if args.haplotype_seqs is None:
            # Strategy 1 - Simulate haplotypes from scratch
            # Simulate master sequence and save record to FASTA file
            haplotype_seq = sim_master(length=args.genome_length,
                                       outdir_haps=outdir_haps, seed=seed)
            if args.use_master:
                # Strategy 1.a - Generate haplotype sequences from a master
                # haplotype
                if args.verbose:
                    print_params(args.mutation_rate, args.deletion_rate,
                                 args.insertion_rate)
                haplotype_seqs = [mutate(
                    haplotype_seq, args.mutation_rate, args.deletion_rate,
                    args.insertion_rate, noFR=args.no_FR,
                    del_len=args.deletion_length, seed=seed + i,
                    verbose=args.verbose) for i in range(num_haplotypes)]
            elif args.tree_like:
                # Strategy 1.b - Sample genotypes from the leaves of a perfect
                # tree
                # Max depth of the perfect tree
                max_depth = np.ceil(np.log2(num_haplotypes))
                root = Tree(haplotype_seq, max_depth)
                # Populate tree
                leaf(root, args.mutation_rate, args.deletion_rate,
                     args.insertion_rate, args.no_FR, args.deletion_length, seed,
                     verbose=args.verbose)
                sequences = root.get_leaves()

                np.random.seed(seed)
                idxs = np.random.randint(0, len(sequences) - 1,
                                         size=num_haplotypes)
                haplotype_seqs = [sequences[idx] for idx in idxs]
            else:
                # Strategy 1.b - Generate num_haplotypes sequences randomly,
                # each of length args.genome_length
                # FIXME: Replace hard coded weights
                weights = (('A', 0.245), ('C', 0.245),
                           ('G', 0.245), ('T', 0.245), ('-', 0.02))
                haplotype_seqs = [sim_haplotypes(
                    length=args.genome_length, weights=weights, s=seed)
                    for _x in range(num_haplotypes)]

            # Save each haplotype sequence to a separate file
            write_fasta(haplotype_seqs, outdir_haps)
        else:
            # Strategy 2 - Simulate haplotypes from a file providing the
            # underlying
            # haplotype(s)
            if args.use_master:
                # Strategy 2.a - From the first sequence provided generate
                # num_haplotypes sequences by mutating the original sequence.
                header, haplotype_seq = read_fasta(args.haplotype_seqs)

                if args.verbose:
                    print_params(args.mutation_rate, args.deletion_rate,
                                 args.insertion_rate)
                haplotype_seqs = [mutate(
                    haplotype_seq[0], args.mutation_rate, args.deletion_rate,
                    args.insertion_rate, noFR=args.no_FR,
                    del_len=args.deletion_length, seed=seed + i,
                    verbose=args.verbose) for i in range(num_haplotypes)]
                write_fasta(haplotype_seqs, outdir_haps)
            elif args.tree_like:
                # Strategy 2.b - Sample genotypes from the leaves of a perfect
                # tree
                # Max depth of the perfect tree
                header, haplotype_seq = read_fasta(args.haplotype_seqs)

                max_depth = np.ceil(np.log2(num_haplotypes))
                root = Tree(haplotype_seq[0], max_depth)
                # Populate tree
                leaf(root, args.mutation_rate, args.deletion_rate,
                     args.insertion_rate, args.no_FR, args.deletion_length,
                     seed, verbose=args.verbose)
                sequences = root.get_leaves()

                np.random.seed(seed)
                idxs = np.random.randint(0, len(sequences) - 1,
                                         size=num_haplotypes)
                haplotype_seqs = [sequences[idx] for idx in idxs]
                write_fasta(haplotype_seqs, outdir_haps)
            else:
                # Strategy 2.c - Split haplotypes into different files
                # (implicitly, it is assumed the number of haplotypes sequences
                # is larger than 1)
                idx = 0
                aux = []
                output_file = ''
                haplotype_id = ''
                with open(args.haplotype_seqs, 'r') as infile:
                    # The following is needed because the number of characters
                    # per line in a FASTA file often do not exceed a certain number
                    for line in infile:
                        record = line.rstrip()
                        if record and record[0] == '>':
                            if idx > 0:
                                with open(output_file, 'w') as outfile:
                                    outfile.write(haplotype_id + '\n')
                                    outfile.write(''.join(aux))
                                aux = []
                            haplotype_id = record
                            output_file = os.path.join(outdir_haps, ''.join(
                                ("haplotype", str(idx), ".fasta")))
                            idx += 1
                        else:
                            aux.append(record)

                output_file = os.path.join(outdir_haps, ''.join(
                    ("haplotype", str(idx - 1), ".fasta")))
                with open(output_file, 'w') as outfile:
                    outfile.writelines([haplotype_id, '\n', ''.join(aux)])
                num_haplotypes = idx
                haplotype_file = os.path.join(outdir_haps, "haplotypes.fasta")
                shutil.copyfile(args.haplotype_seqs, haplotype_file)

    if args.output == "reads" or args.output == "all":
        fragment_mean, fragment_sd = [int(x)
                                      for x in args.fragment_size.split(",")]

        coverage = [int(x) for x in args.coverage.split(",")]
        if len(coverage) > 1:
            assert len(coverage) == num_haplotypes, (
                    "More than one value for coverage specified, but it does not "
                    "coincide with the number of haplotypes")

        if args.freq_dstr == 'unif':
            if len(coverage) > 1:
                print("More than one value for coverage specified, only first "
                      "value is used")

            haplotype_file = os.path.join(outdir_haps, "haplotypes.fasta")
            print(f"Reading file containing sequences of underlying haplotypes "
                  f"from {haplotype_file}")
            tmp_file = os.path.join(outdir_haps, "tmp.fasta")
            shutil.copyfile(haplotype_file, tmp_file)
            # remove true deletions (if present) before generating reads
            sed('s/-//g', tmp_file, verbose=args.verbose)
            outprefix = os.path.join(outdir_reads, "simreads_R")
            sim_reads(args.art, haplotype_seq=tmp_file, coverage=coverage[0],
                      read_len=args.read_length, fragment_mean=fragment_mean,
                      fragment_sd=fragment_sd, outprefix=outprefix,
                      paired=args.paired, highQ=args.quality,
                      num_reads=args.num_reads, seed=seed, verbose=args.verbose)
            os.remove(tmp_file)

            # Make headers compatible with output from Illumina platforms
            # (expected by ngshmmalign)
            if args.paired:
                sed('/^@.*-[0-9]*\/1$/ s/\/1$/ 1:N:0:5/',
                    ''.join((outprefix, '1.fq')), verbose=args.verbose)
                sed('/^@.*-[0-9]*\/2$/ s/\/2$/ 2:N:0:5/',
                    ''.join((outprefix, '2.fq')), verbose=args.verbose)

            # Rename output file
            if args.paired:
                os.rename(''.join((outprefix, '1.fq')),
                          ''.join((outprefix, '1.fastq')))
                os.rename(''.join((outprefix, '2.fq')),
                          ''.join((outprefix, '2.fastq')))
            else:
                os.rename(''.join((outprefix, '.fq')),
                          ''.join((outprefix, '1.fastq')))

        elif args.freq_dstr == 'geom' or args.freq_dstr == 'dirichlet' or args.freq_dstr == 'cust':

            if args.freq_dstr == 'geom':
                if len(coverage) > 1:
                    print(
                        "More than one value for coverage specified, only first "
                        "value is used")
                freqs = [args.geom_ratio**(i + 1) for i in range(num_haplotypes)]
                freqs = np.asarray(freqs)
                freqs = freqs / np.sum(freqs)
                np.set_printoptions(precision=4)
                if args.verbose:
                    print("Relative abundances: ", freqs)
                coverage = freqs * coverage[0]
                coverage = coverage.astype(int)
            elif args.freq_dstr == 'dirichlet':
                if len(coverage) > 1:
                    print(
                        "More than one value for coverage specified, only first "
                        "value is used")
                if args.dirichlet_alpha is None:
                    alpha = np.ones(num_haplotypes)
                else:
                    alpha = [float(x) for x in args.dirichlet_alpha.split(",")]
                    if len(alpha) == 1:
                        alpha = np.repeat(alpha, num_haplotypes)
                    assert len(alpha) == num_haplotypes, (
                            "The number of Dirichlet parameters and number of "
                            "haplotypes does not coincide")
                np.random.seed(seed)
                freqs = np.random.dirichlet(alpha)
                np.set_printoptions(precision=4)
                if args.verbose:
                    print("Relative abundances: ", freqs)
                coverage = freqs * coverage[0]
                coverage = coverage.astype(int)
                # write to output
                fasta_record = collections.namedtuple("fasta_record", "id seq")
                output_file = os.path.join(
                    outdir_haps, "haplotype_frequencies.fasta")
                with open(output_file, 'w') as outfile:
                    for i in range(num_haplotypes):
                        haplotype_id = ''.join(("haplotype", str(i)))
                        line = fasta_record(id=haplotype_id, seq=freqs[i])
                        outfile.write(">{}\n{}\n".format(line.id, line.seq))
            elif args.freq_dstr == 'cust':
                print("Not implemented yet!")
                sys.exit()

            if args.paired:
                outfiles_R1 = []
                outfiles_R2 = []
            else:
                outfiles = []
            for idx in range(num_haplotypes):
                haplotype_file = os.path.join(
                    outdir_haps, ''.join(("haplotype", str(idx), ".fasta")))
                print(f"Reading file containing sequence of haplotype {idx} "
                      f"from {haplotype_file}")
                # remove true deletions (if present) before generating reads
                sed('s/-//g', haplotype_file, verbose=args.verbose)
                outprefix = os.path.join(
                    outdir_reads, ''.join(("reads_hap", str(idx), "_R")))
                if args.paired:
                    outfiles_R1.append(''.join((outprefix, "1.fq")))
                    outfiles_R2.append(''.join((outprefix, "2.fq")))
                else:
                    outfiles.append(''.join((outprefix, ".fq")))

                sim_reads(args.art, haplotype_seq=haplotype_file,
                          coverage=coverage[idx], read_len=args.read_length,
                          fragment_mean=fragment_mean, fragment_sd=fragment_sd,
                          outprefix=outprefix, paired=args.paired,
                          highQ=args.quality, num_reads=args.num_reads,
                          seed=seed, verbose=args.verbose)

            if args.paired:
                sh.cat(outfiles_R1, _out=os.path.join(
                    outdir_reads, "simreads_R1.fastq"))
                sh.cat(outfiles_R2, _out=os.path.join(
                    outdir_reads, "simreads_R2.fastq"))
                for idx in range(len(outfiles_R1)):
                    os.remove(outfiles_R1[idx])
                    os.remove(outfiles_R2[idx])
            else:
                sh.cat(outfiles, _out=os.path.join(
                    outdir_reads, "simreads_R1.fastq"))
                for f in outfiles:
                    os.remove(f)

            # Make headers compatible with output from Illumina platforms
            # (expected by ngshmmalign)
            if args.paired:
                sed('/^@.*-[0-9]*\/1$/ s/\/1$/ 1:N:0:5/',
                    os.path.join(outdir_reads, "simreads_R1.fastq"),
                    verbose=args.verbose)
                sed('/^@.*-[0-9]*\/2$/ s/\/2$/ 2:N:0:5/',
                    os.path.join(outdir_reads, "simreads_R2.fastq"),
                    verbose=args.verbose)

if __name__ == '__main__':
    main()
