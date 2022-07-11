#!/bin/env perl
# extract sequence and info for defined inserts
use strict;
use warnings;

use Getopt::Long;
use Bio::SeqIO;

my $usage = "Usage: $0 -i INSERT-POS.BED -s LONGREAD-SEQ.FASTQ -fq FASTQ-OUTFILE -fa FASTA-OUTFILE -info TSV-OUTFILE";

# get opts
my $insert_file;
my $read_file;
my $fq_outfile;
my $fa_outfile;
my $info_outfile;

@ARGV >= 10 or die $usage;

GetOptions(
"i=s" => \$insert_file,
"s=s" => \$read_file,
"fq=s" => \$fq_outfile,
"fa=s" => \$fa_outfile,
"info=s" => \$info_outfile)
or die "Error in command line arguments, usage: $usage";

defined($insert_file) && defined($read_file) && defined($fq_outfile) && defined($fa_outfile) && defined($info_outfile) || die $usage;

# open input
open(INS, "<$insert_file") || die "Unable to open $insert_file: $!";
open(FQ, "<$read_file") || die "Unable to open $read_file: $!";

open(FQO, ">$fq_outfile") || die "Unable to write to $fq_outfile: $!";
my $fao = new Bio::SeqIO(-file => ">$fa_outfile", -format => 'fasta');
open(INFO, ">$info_outfile") || die "Unable to write to $info_outfile: $!";

$fao->width(100); # set a line-width

my @headers = qw(insert_id insert_chr insert_start insert_end insert_strand insert_len insert_left insert_right insert_detect_type);
print INFO join("\t", @headers), "\n";

# get merged
my %rid2mid;
my %mid2info;
while(my $line = <INS>) {
	chomp $line;
	my ($rid, $from, $to, $inames, $qlen, $strand) = split(/\t/, $line);
	my $len = $to - $from;
	my ($m_chr, $m_strand, $m_start, $m_end);
	foreach my $iname (split(/,/, $inames)) {
		my @iloc = split(/:/, $iname);
		($m_chr, $m_strand) = @iloc[0,1];
		$m_start = $iloc[2] if(!defined $m_start || $iloc[2] < $m_start);
		$m_end = $iloc[3] if(!defined $m_end || $iloc[3] > $m_start);
	}

  my $m_len = $to - $from;
	my $m_left = $from;
	my $m_right = $qlen - $to;
	if($strand eq '-') {
		($m_left, $m_right) = ($m_right, $m_left);
	}

	my $mtype = $m_left > 0 && $m_right > 0 ? 'complete' : 'incomplete';
	my $mop = $m_left > 0 && $m_right > 0 ? 'I' : 'S';
	my $mid = "$rid:$m_chr:$m_strand:$m_start:$m_end:$m_left" . 'L' . ':' . $len . $mop . ':' . $m_right . 'R';

  push(@{$rid2mid{$rid}}, $mid);
	@{$mid2info{$mid}} = ($rid, $from, $to, $inames, $qlen, $strand);
	print INFO "$mid\t$m_chr\t$m_start\t$m_end\t$m_strand\t$m_len\t$m_left\t$m_right\t$mtype\n";
}

# search long read and output
while(my $def = <FQ>) {
	chomp $def;
	my ($rid) = $def =~ /^@(\S+)/;
	my $rseq = <FQ>; chomp $rseq;
	my $rsep = <FQ>; chomp $rsep;
	my $rqual = <FQ>; chomp $rqual;

	foreach my $mid (@{$rid2mid{$rid}}) {
		my ($rid, $from, $to, $inames, $rlen, $strand) = @{$mid2info{$mid}};
		my $mseq = $rseq;
		my $mqual = $rqual;
		if($strand eq '-') { # revcom
			$mseq = reverse($mseq);
			$mqual = reverse($mqual);
			$mseq =~ tr/ACGTUacgtu/TGCAAtgcaa/;
		}
		$mseq = substr($mseq, $from, $to - $from);
		$mqual = substr($mqual, $from, $to - $from);

# output
    print FQO "@", "$mid\n$mseq\n$rsep\n$mqual\n";
		$fao->write_seq(new Bio::Seq(-seq => $mseq, -display_id => $mid, -desc => $inames));
	}
}


close(INS);
close(FQ);
close(FQO);
$fao->close();
close(INFO);
