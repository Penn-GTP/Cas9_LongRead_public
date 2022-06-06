#!/bin/env perl
# This script used to get the sequence, region, and annotation files from the VECTOR GenBank file for remap annotation purpose
use strict;
use warnings;

use Bio::SeqIO;
use Bio::Tools::GFF;
my $gff_formatter = new Bio::Tools::GFF(-gff_version => 3); # GFF formatter

my $usage = "Usage: $0 GenBank-INFILE SEQ-OUTFILE REGION-OUTFILE ANNO-OUTFILE";
my $infile = shift or die $usage;
my $seq_outfile = shift or die $usage;
my $reg_outfile = shift or die $usage;
my $anno_outfile = shift or die $usage;

# open GenBank input
my $in = new Bio::SeqIO(-file => "<$infile", -format => 'genbank');

# open FASTA seq output
my $out = new Bio::SeqIO(-file => ">$seq_outfile", -format => 'fasta');

# open BED output
open(REG, ">$reg_outfile") || die "Unable to write to: $!";

# open GFF output
open(ANNO, ">$anno_outfile") || die "Unable to write to: $!";
# add headers
print ANNO "##gff-version 3\n";
# customized trackline
print ANNO "##displayName=label\n";

# read each GenBank seq
while(my $seqobj = $in->next_seq()) {
	my $seqID = $seqobj->display_id();
	my $len = $seqobj->length();
	print REG "$seqID\t0\t$len\n";

	my $featID = 0;
	$out->write_seq($seqobj);
# scan each FEATURE entity
	foreach my $feat ($seqobj->get_SeqFeatures()) {
# add uniq IDs
		$feat->add_tag_value('ID', "$seqID." . (++$featID));
# add a name for featureCount
    my ($name) = $feat->has_tag('label') ? $feat->get_tag_values('label') : ($feat->primary_tag());
		$feat->add_tag_value('Name', $name);
		print ANNO $feat->gff_string($gff_formatter), "\n";
	}
}

$in->close();
$out->close();
close(ANNO);
