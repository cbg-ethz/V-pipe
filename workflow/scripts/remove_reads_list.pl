#!/usr/bin/env perl 

use strict;
use warnings;

# reads sequence ids from file from first command line argument
# then iterates over stdin and removes sequences for the given ids
# the sequences from stdin are organized in groups of three lines, the first line
# contains the sequence id

my %h; # hash with 'headers' to filter-out as keys

{
	my $listfile = shift @ARGV;
	open my $fh, '<', $listfile
		or die "cannot open ${listfile}";

	# build a hash out of all (unique) seq ID
	%h = map { chomp; ('@' . $_) => 1 } <$fh>;

	close $fh;
}
print STDERR "filtering out: ", scalar keys %h, " seq id\n";

while (<>) { # loop through headers...
	my @s = map { scalar <> } 1 .. 3; # ...fetch the remainer 3 lines of read: bp sequence, separator, phred quality scores
	print $_, @s
		unless(exists $h{(split)[0]}); # consider header's first part ('@' + seq ID), ignore the second half ("comments", e.g.: Illumina indexes)
}
