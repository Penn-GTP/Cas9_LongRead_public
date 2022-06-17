#!/bin/env perl
# extract sequence and info of off-target insertions/soft-clips
use strict;
use warnings;

use Getopt::Long;
use Bio::SeqIO;

my $min_insert = 20;
my $usage = "Usage: $0 -t TARGET-BEDFILE -i INFILE -fq FASTQ-OUTFILE -fa FASTA-OUTFILE -info TSV-OUTFILE [--min-insert $min_insert]";

# get opts
my $target_file;
my $infile;
my $fq_outfile;
my $fa_outfile;
my $info_outfile;

@ARGV >= 10 or die $usage;

GetOptions(
"t=s" => \$target_file,
"i=s" => \$infile,
"fq=s" => \$fq_outfile,
"fa=s" => \$fa_outfile,
"info=s" => \$info_outfile,
"min-insert=i" => \$min_insert)
or die "Error in command line arguments, usage: $usage";

defined($target_file) && defined($infile) && defined($fq_outfile) && defined($fa_outfile) && defined($info_outfile) || die $usage;

# open input
open(BED, "<$target_file") || die "Unable to open $target_file: $!";
open(IN, "samtools view $infile |") || die "Unable to open samtools with samtools: $!";
open(FQO, ">$fq_outfile") || die "Unable to write to $fq_outfile: $!";
my $fao = new Bio::SeqIO(-file => ">$fa_outfile", -format => 'fasta', -alphabet => 'dna');
open(INFO, ">$info_outfile") || die "Unable to write to $info_outfile: $!";

$fao->width(100); # set a line-width

my @headers = qw(insert_id insert_chr insert_pos insert_strand insert_len insert_left insert_right insert_rel_pos insert_detect_type);
print INFO join("\t", @headers), "\n";

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
		elsif($op eq 'I') { # insert
			if(!($chr eq $rname && $start <= $insert_start && $insert_start + 1 <= $end) && $len >= $min_insert) { # current insert location is off target region
				my $insert_left = $insert_from;
				my $insert_right = $qlen - $insert_from - $len;
				my $insert_seq = substr($seq, $insert_from, $len);
				my $insert_qual = substr($qual, $insert_from, $len);
				die if(length($insert_seq) != length($insert_qual));
				if($strand eq '-') { # revcom seq and reverse qual
					$insert_seq = reverse($insert_seq);
					$insert_qual = reverse($insert_qual);
					$insert_seq =~ tr/ACGTUacgtu/TGCAAtgcaa/;
					($insert_left, $insert_right) = ($insert_right, $insert_left); # swap left-right size
				}
				my $insert_id = "$qname:$chr:$start:$end:$target_name:$chr:$strand:$insert_start:$insert_left" . 'L:' . "$len$op:$insert_right" . 'R';
				my $rel_pos = $chr eq $rname ? $insert_start - $start : 'inf';
				my $detect_type = $op eq 'I' ? 'complete' : 'incomplete';
				print FQO "\@$insert_id\n$insert_seq\n+\n$insert_qual\n";
				$fao->write_seq(new Bio::Seq(-display_id => $insert_id, -seq => $insert_seq));
				print INFO "$insert_id\t$rname\t$insert_start\t$strand\t$len\t$insert_left\t$insert_right\t$rel_pos\t$detect_type\n";
			}
			$insert_from += $len; # update
		}
		elsif($op eq 'S') { # soft-clip
		  $insert_from += $len;
		}
		else {
			next;
		}
	}
}


close(BED);
close(IN);
close(FQO);
$fao->close();
