#!/bin/env perl
# extract sequence and info of on-target insertions/soft-clips
use strict;
use warnings;

use Getopt::Long;

my $min_insert = 20;
my $max_dist = 500;
my $usage = "Usage: $0 -t TARGET-BEDFILE -i INFILE -o OUTFILE [--min-insert $min_insert] [--max-dist $max_dist]";

# get opts
my $target_file;
my $infile;
my $outfile;

@ARGV >= 6 or die $usage;

GetOptions(
"t=s" => \$target_file,
"i=s" => \$infile,
"o=s" => \$outfile,
"min-insert=i" => \$min_insert,
"max-dist=i" => \$max_dist)
or die "Error in command line arguments, usage: $usage";

defined($target_file) && defined($infile) && defined($outfile) || die $usage;

# open input
open(BED, "<$target_file") || die "Unable to open $target_file: $!";
open(IN, "samtools view $infile |") || die "Unable to open samtools with samtools: $!";
open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";

# read in target site
my $loc = <BED>;
chomp $loc;
my ($chr, $start, $end, $target_name) = split(/\t/, $loc);

while(my $line = <IN>) {
	chomp $line;
	my ($qname, $flag, $rname, $pos, $mapQ, $cigar, $rnext, $pnext, $tlen, $seq, $qual, @tags) = split(/\t/, $line);
	my $qlen = length $seq;
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
		elsif($op eq 'D') {
			$insert_start += $len;
		}
		elsif($op eq 'S' || $op eq 'I') { # insert or 5'/3' end clip region
			if($chr eq $rname && $start - $max_dist <= $insert_start && $insert_start + 1 <= $end + $max_dist && $len >= $min_insert) { # current insert location is on target region
				my $insert_to = $insert_from + $len;
				my $insert_end = $insert_start + 1;
				my $insert_id = "$rname:$strand:$insert_start:$insert_end:$len$op";
				print OUT "$qname\t$insert_from\t$insert_to\t$insert_id\t$qlen\t$strand\n";
			}
			$insert_from += $len; # update
		}
		else {
			next;
		}
	}
}


close(IN);
close(OUT);
