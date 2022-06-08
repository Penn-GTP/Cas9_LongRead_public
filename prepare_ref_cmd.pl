#!/bin/env perl
# Prepare sh script building Cas9LongRead REF
our $VERSION = v1.1;
our $ENV_FILE = 'set_ref_env.sh';

use strict;
use warnings;
use lib '/project/gtplab/pipeline/Cas9_LongRead';
use Cas9LongReadExpDesign;

my $usage = "Usage: perl $0 DESIGN-FILE BASH-OUTFILE";
my $sh_path = '/bin/bash';
my $vec_anno_script = 'get_vector_info.pl';
my $samtools = 'samtools';

my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $design = new Cas9LongReadExpDesign($infile);
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
	print STDERR "Error: VECTOR_DIR $VEC_DIR not exists\n";
	exit;
}

if(!(-e $SCRIPT_DIR && -d $SCRIPT_DIR)) {
	print STDERR "Error: SCRIPT_DIR $SCRIPT_DIR not exists\n";
	exit;
}

if(!(-e $WORK_DIR)) {
	print STDERR "Warning: WORK_DIR $WORK_DIR not exists, creating now\n";
  mkdir($WORK_DIR, 0750) || die "Unable to mkdir $WORK_DIR: $!\n";
}

# check global options
if(!($NGS_ALIGNER eq 'minimap2')) {
	print STDERR "Error: only support 'minimap2' for the NGS_ALIGNER option\n";
	exit;
}

open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";
# write header
print OUT "#!$sh_path\n";
# set env
print OUT "source $SCRIPT_DIR/$ENV_FILE\n\n";

foreach my $sample ($design->get_sample_names()) {
# prepare get nuclease vec info cmd
	if($design->sample_opt($sample, 'nuclease_gb')) {
		my $in = $design->sample_opt($sample, 'nuclease_gb');
		my $seq_out = $design->get_sample_nuclease_vec_seq($sample);
		my $reg_out = $design->get_sample_nuclease_vec_region($sample);
		my $anno_out = $design->get_sample_nuclease_vec_anno($sample);
		my $cmd = "$SCRIPT_DIR/$vec_anno_script $VEC_DIR/$in $VEC_DIR/$seq_out $VEC_DIR/$reg_out $VEC_DIR/$anno_out";
		if(!(-e "$VEC_DIR/$seq_out" && -e "$VEC_DIR/$reg_out" && -e "$VEC_DIR/$anno_out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: nuclease vec seq and annotation file exist, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare get donor vec info cmd
  if($design->sample_opt($sample, 'donor_gb')) {
		my $in = $design->sample_opt($sample, 'donor_gb');
		my $seq_out = $design->get_sample_donor_vec_seq($sample);
		my $reg_out = $design->get_sample_donor_vec_region($sample);
		my $anno_out = $design->get_sample_donor_vec_anno($sample);
		my $cmd = "$SCRIPT_DIR/$vec_anno_script $VEC_DIR/$in $VEC_DIR/$seq_out $VEC_DIR/$reg_out $VEC_DIR/$anno_out";
		if(!(-e "$VEC_DIR/$seq_out" && -e "$VEC_DIR/$reg_out" && -e "$VEC_DIR/$anno_out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: donor vec seq and annotation file exist, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare get trans vec info cmd
  if($design->sample_opt($sample, 'trans_gb')) {
		my $in = $design->sample_opt($sample, 'trans_gb');
		my $seq_out = $design->get_sample_trans_vec_seq($sample);
		my $reg_out = $design->get_sample_trans_vec_region($sample);
		my $anno_out = $design->get_sample_trans_vec_anno($sample);
		my $cmd = "$SCRIPT_DIR/$vec_anno_script $VEC_DIR/$in $VEC_DIR/$seq_out $VEC_DIR/$reg_out $VEC_DIR/$anno_out";
		if(!(-e "$VEC_DIR/$seq_out" && -e "$VEC_DIR/$reg_out" && -e "$VEC_DIR/$anno_out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: donor vec seq and annotation file exist, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare get helper vec info cmd
  if($design->sample_opt($sample, 'helper_gb')) {
		my $in = $design->sample_opt($sample, 'helper_gb');
		my $seq_out = $design->get_sample_helper_vec_seq($sample);
		my $reg_out = $design->get_sample_helper_vec_region($sample);
		my $anno_out = $design->get_sample_helper_vec_anno($sample);
		my $cmd = "$SCRIPT_DIR/$vec_anno_script $VEC_DIR/$in $VEC_DIR/$seq_out $VEC_DIR/$reg_out $VEC_DIR/$anno_out";
		if(!(-e "$VEC_DIR/$seq_out" && -e "$VEC_DIR/$reg_out" && -e "$VEC_DIR/$anno_out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: donor vec seq and annotation file exist, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare combine vec seq and index cmd
	{
    my $out = $design->get_sample_vec_seq($sample);
    my @infiles;
		if($design->sample_opt($sample, 'nuclease_gb')) {
			push(@infiles, "$VEC_DIR/" . $design->get_sample_nuclease_vec_seq($sample));
		}
		if($design->sample_opt($sample, 'donor_gb')) {
			push(@infiles, "$VEC_DIR/" . $design->get_sample_donor_vec_seq($sample));
		}
		if($design->sample_opt($sample, 'trans_gb')) {
			push(@infiles, "$VEC_DIR/" . $design->get_sample_trans_vec_seq($sample));
		}
		if($design->sample_opt($sample, 'helper_gb')) {
			push(@infiles, "$VEC_DIR/" . $design->get_sample_helper_vec_seq($sample));
		}

		my $cmd = "cat " . join(" ", @infiles) . " > $VEC_DIR/$out";

# add faidx cmd
		$cmd .= "\n$samtools faidx $VEC_DIR/$out";

		if(!-e "$VEC_DIR/$out") {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $VEC_DIR/$out exists, won't override\n";
			$cmd =~ s/\n/\n# /sg; # comment out to intermediate lines
			print OUT "# $cmd\n";
		}
	}

# prepare combine vec anno cmd
	{
    my $out = $design->get_sample_vec_anno($sample);
    my @infiles;
		if($design->sample_opt($sample, 'nuclease_gb')) {
			push(@infiles, "$VEC_DIR/" . $design->get_sample_nuclease_vec_anno($sample));
		}
		if($design->sample_opt($sample, 'donor_gb')) {
			push(@infiles, "$VEC_DIR/" . $design->get_sample_donor_vec_anno($sample));
		}
		if($design->sample_opt($sample, 'trans_gb')) {
			push(@infiles, "$VEC_DIR/" . $design->get_sample_trans_vec_anno($sample));
		}
		if($design->sample_opt($sample, 'helper_gb')) {
			push(@infiles, "$VEC_DIR/" . $design->get_sample_helper_vec_anno($sample));
		}

		my $cmd = "cat " . join(" ", @infiles) . " > $VEC_DIR/$out";

		if(!-e "$VEC_DIR/$out") {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $VEC_DIR/$out exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

	print OUT "\n";
}

close(OUT);
# change to exacutable
chmod 0750, $outfile;
