#!/bin/env perl
# Prepare bash script for generating per-sample processing stats
our $VERSION = v1.1;

use strict;
use warnings;
use lib '/project/gtplab/pipeline/Cas9_LongRead';
use Cas9OntSeqExpDesign;

my $usage = "Usage: perl $0 DESIGN-FILE BASH-OUTFILE";
#my $sh_path = '/bin/bash';
my $samtools = 'samtools';
my @headers = qw(sample_name total_read ref_mapped ref_enrich ref_target ref_insert insert_mapped insert_nuclease_mapped insert_donor_mapped insert_trans_mapped insert_helper_mapped insert_ref2_mapped insert_vec2_mapped);

my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $design = new Cas9OntSeqExpDesign($infile);
my $NUM_PROC = $design->get_global_opt('NUM_PROC');
my $BASE_DIR = $design->get_global_opt('BASE_DIR');
my $SCRIPT_DIR = $design->get_global_opt('SCRIPT_DIR');
my $WORK_DIR = $design->get_global_opt('WORK_DIR');
my $VEC_DIR = $design->get_global_opt('VEC_DIR');

# check required directories
if(!(-e $BASE_DIR && -d $BASE_DIR)) {
  print STDERR "Error: $BASE_DIR not exists\n";
  exit;
}

if(!(-e $VEC_DIR && -d $VEC_DIR)) {
  print STDERR "Error: $VEC_DIR not exists\n";
  exit;
}

if(!(-e $SCRIPT_DIR && -d $SCRIPT_DIR)) {
  print STDERR "Error: $SCRIPT_DIR not exists\n";
  exit;
}

if(!(-e $WORK_DIR)) {
  print STDERR "Error: WORK_DIR $WORK_DIR not exists\n";
  exit;
}

open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";
# write header
print OUT join("\t", @headers), "\n";

foreach my $sample ($design->get_sample_names()) {
	print STDERR "gathering stats for $sample\n";
# get total read
  my $total_read;
	{
		my $in = $design->sample_opt($sample, 'read_fastq');
		$total_read = $in =~ /\.gz$/ ? `zcat $in | wc -l` : `cat $in | wc -l`;
		chomp $total_read;
		$total_read /= 4;
	}

	my $min_mapQ = $design->sample_opt($sample, 'min_mapQ');
# get ref mapped
  my $ref_mapped;
	{
		my $in = $design->get_sample_ref_map_file($sample);
		$ref_mapped = `$samtools view -F 0x4 -q $min_mapQ $WORK_DIR/$in | cut -f1 | sort -u | wc -l`;
		chomp $ref_mapped;
	}

# get ref enrich
  my $ref_enrich;
	{
		my $in = $design->get_sample_ref_map_enrich_sort_file($sample);
		$ref_enrich = `$samtools view $BASE_DIR/$in | cut -f1 | sort -u | wc -l`;
		chomp $ref_enrich;
	}

# get ref target
  my $ref_target;
	{
		my $in = $design->get_sample_ref_map_target_sort_file($sample);
		$ref_target = `$samtools view $BASE_DIR/$in | cut -f1 | sort -u | wc -l`;
		chomp $ref_target;
	}

# get ref insert
  my $ref_insert;
	{
		my $in = $design->get_sample_target_insert_fastq($sample);
		$ref_insert = $in =~ /\.gz$/ ? `zcat $BASE_DIR/$in | wc -l` : `cat $BASE_DIR/$in | wc -l`;
		chomp $ref_insert;
		$ref_insert /= 4;
	}

# get vec mapped
	my $vec_mapped = 0;
	{
		my $in = $design->get_sample_vec_sorted_file($sample);
		$vec_mapped = `$samtools view | cut -f1 | sort -u | wc -l`;
		chomp $vec_mapped;
	}

# get nuclease vec mapped
	my $nuclease_mapped = 0;
	if($design->sample_opt($sample, 'nuclease_gb')) {
		my $in = $design->get_sample_vec_sorted_file($sample);
		my $bed = $design->get_sample_nuclease_vec_region($sample);
		$nuclease_mapped = `$samtools view -L $VEC_DIR/$bed $BASE_DIR/$in | cut -f1 | sort -u | wc -l`;
		chomp $nuclease_mapped;
	}

# get donor vec mapped
	my $donor_mapped = 0;
	if($design->sample_opt($sample, 'donor_gb')) {
		my $in = $design->get_sample_vec_sorted_file($sample);
		my $bed = $design->get_sample_donor_vec_region($sample);
		$donor_mapped = `$samtools view -L $VEC_DIR/$bed $BASE_DIR/$in | cut -f1 | sort -u | wc -l`;
		chomp $donor_mapped;
	}

# get trans vec mapped
	my $trans_mapped = 0;
	if($design->sample_opt($sample, 'trans_gb')) {
		my $in = $design->get_sample_vec_sorted_file($sample);
		my $bed = $design->get_sample_trans_vec_region($sample);
		$nuclease_mapped = `$samtools view -L $VEC_DIR/$bed $BASE_DIR/$in | cut -f1 | sort -u | wc -l`;
		chomp $nuclease_mapped;
	}

# get helper vec mapped
	my $helper_mapped = 0;
  if($design->sample_opt($sample, 'helper_gb')) {
		my $in = $design->get_sample_vec_sorted_file($sample);
		my $bed = $design->get_sample_helper_vec_region($sample);
		$donor_mapped = `$samtools view -L $VEC_DIR/$bed $BASE_DIR/$in | cut -f1 | sort -u | wc -l`;
		chomp $donor_mapped;
	}

# get ref2 mapped
	my $ref2_mapped = 'NA';
	if($design->get_sample_opts($sample, 'ref2_db')) {
		my $in = $design->get_sample_ref2_sorted_file($sample);
		$ref2_mapped = `$samtools view $BASE_DIR/$in | cut -f1 | sort -u | wc -l`;
		chomp $ref2_mapped;
	}

# get vec2 mapped
	my $vec2_mapped = 'NA';
	if($design->get_sample_opts($sample, 'vec2_db')) {
		my $in = $design->get_sample_vec2_sorted_file($sample);
		$vec2_mapped = `$samtools view $BASE_DIR/$in | cut -f1 | sort -u | wc -l`;
		chomp $vec2_mapped;
	}

# output
  print OUT "$sample\t$total_read\t$ref_mapped\t$ref_enrich\t$ref_target\t$ref_insert\t$vec_mapped\t$nuclease_mapped\t$donor_mapped\t$trans_mapped\t$helper_mapped\t$ref2_mapped\t$vec2_mapped\n";
}

close(OUT);
