#!/bin/env perl
# Prepare bash script for annotate insert based on vec and optionally rev2 and vec2 mapping
our $VERSION = v1.1;
our $ENV_FILE = 'set_annotate_env.sh';

use strict;
use warnings;
use lib '/project/gtplab/pipeline/Cas9_LongRead';
use Cas9OntSeqExpDesign;

my $usage = "Usage: perl $0 DESIGN-FILE BASH-OUTFILE";
my $sh_path = '/bin/bash';
my $insert_anno_script = 'get_insert_anno.pl';
my $samtools = 'samtools';
my $bedtools = 'bedtools';
my $featureCounts = 'featureCounts';

my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $design = new Cas9OntSeqExpDesign($infile);
my $NUM_PROC = $design->get_global_opt('NUM_PROC');
my $BASE_DIR = $design->get_global_opt('BASE_DIR');
my $SCRIPT_DIR = $design->get_global_opt('SCRIPT_DIR');
my $VEC_DIR = $design->get_global_opt('VEC_DIR');
my $WORK_DIR = $design->get_global_opt('WORK_DIR');
my $NGS_ALIGNER = $design->get_global_opt('NGS_ALIGNER');

# check required directories
if(!(-e $BASE_DIR && -d $BASE_DIR)) {
	print STDERR "Error: BASE_DIR $BASE_DIR not exists\n";
	exit;
}

if(!(-e $VEC_DIR && -d $VEC_DIR)) {
	print STDERR "Error: BASE_DIR $VEC_DIR not exists\n";
	exit;
}

if(!(-e $SCRIPT_DIR && -d $SCRIPT_DIR)) {
	print STDERR "Error: SCRIPT_DIR $SCRIPT_DIR not exists\n";
	exit;
}

if(!(-e $WORK_DIR)) {
	print STDERR "Error: WORK_DIR $WORK_DIR not exists\n";
	exit;
}

if(!($NGS_ALIGNER eq 'minimap2')) {
	print STDERR "Error: only support 'minimap2' NGS aligner\n";
	exit;
}

open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";
# write header
print OUT "#!$sh_path\n";
# set env
print OUT "source $SCRIPT_DIR/$ENV_FILE\n\n";

foreach my $sample ($design->get_sample_names()) {
# prepare format vec cmd
	{
		my $in = $design->get_sample_vec_sorted_file($sample);
		my $out = $design->get_sample_vec_sorted_bed($sample);
		
		my $cmd = "$bedtools bamtobed -i $BASE_DIR/$in -cigar > $BASE_DIR/$out";
		if(!(-e "$BASE_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $BASE_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare format ref2 cmd
  if($design->sample_opt($sample, 'ref2_db')) {
		my $in = $design->get_sample_ref2_sorted_file($sample);
		my $out = $design->get_sample_ref2_sorted_bed($sample);
		
		my $cmd = "$bedtools bamtobed -i $BASE_DIR/$in -cigar > $BASE_DIR/$out";
		if(!(-e "$BASE_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $BASE_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare format vec2 cmd
  if($design->sample_opt($sample, 'vec2_db')) {
		my $in = $design->get_sample_vec2_sorted_file($sample);
		my $out = $design->get_sample_vec2_sorted_bed($sample);
		
		my $cmd = "$bedtools bamtobed -i $BASE_DIR/$in -cigar > $BASE_DIR/$out";
		if(!(-e "$BASE_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $BASE_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare annotate nuclease map cmd
  if($design->sample_opt($sample, 'nuclease_gb')) {
		my $gff = $design->get_sample_nuclease_vec_anno($sample);
		my $in = $design->get_sample_vec_sorted_bed($sample);
		my $out = $design->get_sample_insert_nuclease_vec_anno($sample);
		my $opts = $design->sample_opt($sample, 'vec_anno_opts');

		my $cmd = "$SCRIPT_DIR/$insert_anno_script $VEC_DIR/$gff $BASE_DIR/$in $BASE_DIR/$out $opts";

		if(!(-e "$BASE_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $BASE_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare annotate donor map cmd
  if($design->sample_opt($sample, 'donor_gb')) {
		my $gff = $design->get_sample_donor_vec_anno($sample);
		my $in = $design->get_sample_vec_sorted_bed($sample);
		my $out = $design->get_sample_insert_donor_vec_anno($sample);
		my $opts = $design->sample_opt($sample, 'vec_anno_opts');

		my $cmd = "$SCRIPT_DIR/$insert_anno_script $VEC_DIR/$gff $BASE_DIR/$in $BASE_DIR/$out $opts";

		if(!(-e "$BASE_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $BASE_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare annotate trans map cmd
  if($design->sample_opt($sample, 'trans_gb')) {
		my $gff = $design->get_sample_trans_vec_anno($sample);
		my $in = $design->get_sample_vec_sorted_bed($sample);
		my $out = $design->get_sample_insert_trans_vec_anno($sample);
		my $opts = $design->sample_opt($sample, 'vec_anno_opts');

		my $cmd = "$SCRIPT_DIR/$insert_anno_script $VEC_DIR/$gff $BASE_DIR/$in $BASE_DIR/$out $opts";

		if(!(-e "$BASE_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $BASE_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare annotate helper map cmd
  if($design->sample_opt($sample, 'helper_gb')) {
		my $gff = $design->get_sample_helper_vec_anno($sample);
		my $in = $design->get_sample_vec_sorted_bed($sample);
		my $out = $design->get_sample_insert_helper_vec_anno($sample);
		my $opts = $design->sample_opt($sample, 'vec_anno_opts');

		my $cmd = "$SCRIPT_DIR/$insert_anno_script $VEC_DIR/$gff $BASE_DIR/$in $BASE_DIR/$out $opts";

		if(!(-e "$BASE_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $BASE_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare annotate ref2 map cmd
  if($design->sample_opt($sample, 'ref2_db')) {
		my $gff = $design->sample_opt($sample, 'ref2_gff');
		my $in = $design->get_sample_ref2_sorted_bed($sample);
		my $out = $design->get_sample_insert_ref2_anno($sample);
		my $opts = $design->sample_opt($sample, 'ref2_anno_opts');

		my $cmd = "$SCRIPT_DIR/$insert_anno_script $gff $BASE_DIR/$in $BASE_DIR/$out $opts";

		if(!(-e "$BASE_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $BASE_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare annotate vec2 map cmd
  if($design->sample_opt($sample, 'vec2_db')) {
		my $gff = $design->sample_opt($sample, 'vec2_gff');
		my $in = $design->get_sample_vec2_sorted_bed($sample);
		my $out = $design->get_sample_insert_vec2_anno($sample);
		my $opts = $design->sample_opt($sample, 'vec2_anno_opts');

		my $cmd = "$SCRIPT_DIR/$insert_anno_script $gff $BASE_DIR/$in $BASE_DIR/$out $opts";

		if(!(-e "$BASE_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $BASE_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare count nuclease map cmd
	{
		my $gff = $design->get_sample_nuclease_vec_anno($sample);
		my $in = $design->get_sample_vec_sorted_file($sample);
		my $out = $design->get_sample_insert_nuclease_vec_count($sample);
		open(GFF, "<$VEC_DIR/$gff") || die "Unable to open $VEC_DIR/$gff: $!";
		my @featTypes = get_uniq_feat_types(\*GFF);
		close(GFF);
		my $opts = $design->sample_opt($sample, 'nuclease_count_opts');

		my $cmd = "$featureCounts -a $VEC_DIR/$gff -o $BASE_DIR/$out -M -O -f -t \""
		. join(",", @featTypes) . "\" -T $NUM_PROC $opts $BASE_DIR/$in";

		if(!(-e "$BASE_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $BASE_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare count donor map cmd
	{
		my $gff = $design->get_sample_donor_vec_anno($sample);
		my $in = $design->get_sample_vec_sorted_file($sample);
		my $out = $design->get_sample_insert_donor_vec_count($sample);
		open(GFF, "<$VEC_DIR/$gff") || die "Unable to open $VEC_DIR/$gff: $!";
		my @featTypes = get_uniq_feat_types(\*GFF);
		close(GFF);
		my $opts = $design->sample_opt($sample, 'donor_count_opts');

		my $cmd = "$featureCounts -a $VEC_DIR/$gff -o $BASE_DIR/$out -M -O -f -t \""
		. join(",", @featTypes) . "\" -T $NUM_PROC $opts $BASE_DIR/$in";

		if(!(-e "$BASE_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $BASE_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

  print OUT "\n";
}

close(OUT);
# change to exacutable
chmod 0750, $outfile;

sub get_uniq_feat_types {
	my $fh = shift;
	my %type_seen;
	while(my $line = <$fh>) {
		chomp $line;
		next if($line =~ /^#/);
		my ($type) = (split(/\t/, $line))[2];
		if($type ne 'source') {
			$type_seen{$type}++;
		}
	}
	return keys %type_seen;
}
