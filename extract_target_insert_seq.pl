#!/bin/env perl
# This script used to get the sequence and annotation files from the VECTOR GenBank file for remap annotation purpose
use strict;
use warnings;

use Getopt::Long;

my $min_insert = 44;
my $usage = "Usage: $0 -t TARGET-BEDFILE -i INFILE -o OUTFILE [--min-insert $min_insert]";

# get opts
my $target_file;
my $infile;
my $outfile;

GetOptions(
"t=s" => \$target_file,
"i=s" => \$infile,
"o=s" => \$outfile,
"min-insert=i" => \$min_insert)
or die "Error in command line arguments, usage: $usage";

# open input
open(BED, "<$target_file") || die "Unable to open $target_file: $!";
open(IN, "samtools view $infile |") || die "Unable to open samtools with samtools: $!";
open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";

# read in target site
my $loc = <BED>;
chomp $loc;
my ($chr, $start, $end) = split(/\t/, $loc);

while(my $line = <IN>) {
	chomp $line;
	my ($qname, $flag, $rname, $pos, $mapQ, $cigar, $rnext, $pnext, $tlen, $seq, $qual, @tags) = split(/\t/, $line);
	if($rname ne $chr) {
		print STDERR "Warning: ROI alignment chromosome does not match the target one, ignore\n";
		next;
	}
	my $strand = ($flag & 0x10) ? '-' : '+';
	my $insert_start = $pos - 1; # 0-based
	my $insert_from = 0; # 0-based
	while($cigar =~ /(\d+)([MIDNSHPX=])/g) {
		my ($len, $op) = ($1, $2);
		if($op eq 'M' || $op eq '=' || $op eq 'X') {
			$insert_start += $len;
			$insert_from += $len;
		}
		elsif($op eq 'H' || $op eq 'P') {
			next;
		}
		elsif($op eq 'N') {
			$insert_start += $len;
		}
		elsif($op eq 'S') {
			$insert_from += $len;
		}
		elsif($op eq 'D') {
			$insert_start += $len;
		}
		elsif($op eq 'I') { # insert region
			if($start <= $insert_start && $insert_start + 1 <= $end && $len >= $min_insert) { # current insert location is in the target region and large enough
				my $insert_id = "$qname:$insert_start:$strand:$len$op";
				my $insert_seq = substr($seq, $insert_from, $len);
				my $insert_qual = substr($qual, $insert_from, $len);
				die if(length($insert_seq) != length($insert_qual));
				if($strand eq '-') { # revcom seq and reverse qual
					$insert_seq = reverse($insert_seq);
					$insert_qual = reverse($insert_qual);
					$insert_seq =~ tr/ACGTUacgtu/TGCAAtgcaa/;
				}
				print OUT "\@$insert_id\n$insert_seq\n+\n$insert_qual\n";
			}
			$insert_from += $len; # update
		}
	}
}


close(BED);
close(IN);
close(OUT);
