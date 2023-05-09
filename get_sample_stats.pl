#!/bin/env perl
# Prepare bash script for generating per-sample processing stats
our $VERSION = 'v1.2.1';

use strict;
use warnings;
use lib '/project/gtplab/pipeline/Cas9_LongRead';
use Cas9LongReadExpDesign;

my $usage = "Usage: perl $0 DESIGN-FILE OUTFILE";
#my $sh_path = '/bin/bash';
my $samtools = 'samtools';
my $cmd = "$0 " . join(" ", @ARGV);
my $comments = qq(Sample name\tTotal reads\tHost-mapped reads\tHost-mapped Cas9 enriched region mapped reads\tHost-mapped ARCUS target mapped reads\t) .
qq(On-target inserts\tOn-target complete inserts\tOn-target incomplete inserts\tOn-target inserts vector mapped\tOn-target inserts nuclease-vector mapped\tOn-target inserts donor-vector mapped\tOn-target inserts trans-plasmid mapped\tOn-target inserts helper-plasmid mapped\tOn-target inserts host2 mapped\tOn-target inserts vector2 mapped\tOn-target inserts functional count\tOn-target inserts functional clone\tOn-target inserts functional clone frequency\tOn-target inserts type frequency\tOn-target inserts functional type frequency\t) .
qq(Of-target inserts\tOff-target complete inserts\tOff-target incomplete inserts\tOff-target inserts vector mapped\tOff-target inserts nuclease-vector mapped\tOff-target inserts donor-vector mapped\tOff-target inserts trans-plasmid mapped\tOff-target inserts helper-plasmid mapped\tOff-target inserts host2 mapped\tOff-target inserts vector2 mapped\tOff-target inserts functional count\tOff-target inserts functional clone\tOff-target inserts functional clone frequency\tOff-target inserts type frequency\tOff-target inserts functional type frequency);

my @headers = qw(sample_name total_read ref_mapped ref_enrich ref_target
target_insert target_insert_complete target_insert_incomplete target_insert_vec_mapped target_insert_nuclease_mapped target_insert_donor_mapped target_insert_trans_mapped target_insert_helper_mapped target_insert_ref2_mapped target_insert_vec2_mapped target_insert_functional_count target_insert_functional_clone target_insert_functional_freq
target_insert_type_freq target_insert_functional_type_freq
off_insert off_insert_complete off_insert_incomplete off_insert_vec_mapped off_insert_nuclease_mapped off_insert_donor_mapped off_insert_trans_mapped off_insert_helper_mapped off_insert_ref2_mapped off_insert_vec2_mapped off_insert_functional_count off_insert_functional_clone off_insert_functional_freq
off_insert_type_freq off_insert_functional_type_freq);

my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $design = new Cas9LongReadExpDesign($infile);
my $NUM_PROC = $design->get_global_opt('NUM_PROC');
my $BASE_DIR = $design->get_global_opt('BASE_DIR');
my $SCRIPT_DIR = $design->get_global_opt('SCRIPT_DIR');
my $WORK_DIR = $design->get_global_opt('WORK_DIR');
my $VEC_DIR = $design->get_global_opt('VEC_DIR');
my $FEATURE_TAG = $design->get_global_opt('FEATURE_TAG');

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
print OUT qq(# CMD:"$cmd"\n# VER:$VERSION\n);
print OUT "# $comments\n";
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
		my $in = $design->get_sample_ref_map_enrich_sorted_file($sample);
		$ref_enrich = `$samtools view $BASE_DIR/$in | cut -f1 | sort -u | wc -l`;
		chomp $ref_enrich;
	}

# get ref target
  my $ref_target;
	{
		my $in = $design->get_sample_ref_map_target_sorted_file($sample);
		$ref_target = `$samtools view $BASE_DIR/$in | cut -f1 | sort -u | wc -l`;
		chomp $ref_target;
	}

# get taget insert info
  my ($target_insert, $target_complete, $target_incomplete) = (0, 0, 0);
	{
		my $in = $design->get_sample_target_insert_info($sample);
		open(INFO, "<$BASE_DIR/$in") || die "Unable to open $BASE_DIR/$in: $!";
		<INFO>; # discard header
		while(my $line = <INFO>) {
			chomp $line;
			$target_insert++;
			my @fields = split(/\t/, $line);
			my $detect_type = pop @fields;
			$detect_type eq 'complete' ? $target_complete++ : $target_incomplete++;
		}
	}

# get taget insert info
  my ($off_insert, $off_complete, $off_incomplete) = (0, 0, 0);
	{
		my $in = $design->get_sample_off_insert_info($sample);
		open(INFO, "<$BASE_DIR/$in") || die "Unable to open $BASE_DIR/$in: $!";
		<INFO>; # discard header
		while(my $line = <INFO>) {
			chomp $line;
			$off_insert++;
			my @fields = split(/\t/, $line);
			my $detect_type = pop @fields;
			$detect_type eq 'complete' ? $off_complete++ : $off_incomplete++;
		}
	}

# get vec mapped
	my ($target_vec_mapped, $off_vec_mapped) = (0, 0);
	{
		my $in_target = $design->get_sample_target_insert_vec_sorted_file($sample);
		my $in_off = $design->get_sample_off_insert_vec_sorted_file($sample);
		$target_vec_mapped = `$samtools view $BASE_DIR/$in_target | cut -f1 | sort -u | wc -l`;
		chomp $target_vec_mapped;
		$off_vec_mapped = `$samtools view $BASE_DIR/$in_off | cut -f1 | sort -u | wc -l`;
		chomp $off_vec_mapped;
	}

# get nuclease vec mapped
	my ($target_nuclease_mapped, $off_nuclease_mapped) = (0, 0);
	if($design->sample_opt($sample, 'nuclease_gb')) {
		my $in_target = $design->get_sample_target_insert_vec_sorted_file($sample);
		my $in_off = $design->get_sample_off_insert_vec_sorted_file($sample);
		my $bed = $design->get_sample_nuclease_vec_region($sample);
		$target_nuclease_mapped = `$samtools view -L $VEC_DIR/$bed $BASE_DIR/$in_target | cut -f1 | sort -u | wc -l`;
		$off_nuclease_mapped = `$samtools view -L $VEC_DIR/$bed $BASE_DIR/$in_off | cut -f1 | sort -u | wc -l`;
		chomp $target_nuclease_mapped;
		chomp $off_nuclease_mapped;
	}

# get donor vec mapped
	my ($target_donor_mapped, $off_donor_mapped) = (0, 0);
	if($design->sample_opt($sample, 'donor_gb')) {
		my $in_target = $design->get_sample_target_insert_vec_sorted_file($sample);
		my $in_off = $design->get_sample_off_insert_vec_sorted_file($sample);
		my $bed = $design->get_sample_donor_vec_region($sample);
		$target_donor_mapped = `$samtools view -L $VEC_DIR/$bed $BASE_DIR/$in_target | cut -f1 | sort -u | wc -l`;
		$off_donor_mapped = `$samtools view -L $VEC_DIR/$bed $BASE_DIR/$in_off | cut -f1 | sort -u | wc -l`;
		chomp $target_donor_mapped;
		chomp $off_donor_mapped;
	}

# get trans vec mapped
	my ($target_trans_mapped, $off_trans_mapped) = (0, 0);
	if($design->sample_opt($sample, 'trans_gb')) {
		my $in_target = $design->get_sample_target_insert_vec_sorted_file($sample);
		my $in_off = $design->get_sample_off_insert_vec_sorted_file($sample);
		my $bed = $design->get_sample_trans_vec_region($sample);
		$target_trans_mapped = `$samtools view -L $VEC_DIR/$bed $BASE_DIR/$in_target | cut -f1 | sort -u | wc -l`;
		$off_trans_mapped = `$samtools view -L $VEC_DIR/$bed $BASE_DIR/$in_off | cut -f1 | sort -u | wc -l`;
		chomp $target_trans_mapped;
		chomp $off_trans_mapped;
	}

# get helper vec mapped
	my ($target_helper_mapped, $off_helper_mapped) = (0, 0);
	if($design->sample_opt($sample, 'helper_gb')) {
		my $in_target = $design->get_sample_target_insert_vec_sorted_file($sample);
		my $in_off = $design->get_sample_off_insert_vec_sorted_file($sample);
		my $bed = $design->get_sample_helper_vec_region($sample);
		$target_helper_mapped = `$samtools view -L $VEC_DIR/$bed $BASE_DIR/$in_target | cut -f1 | sort -u | wc -l`;
		$off_helper_mapped = `$samtools view -L $VEC_DIR/$bed $BASE_DIR/$in_off | cut -f1 | sort -u | wc -l`;
		chomp $target_helper_mapped;
		chomp $off_helper_mapped;
	}

# get ref2 mapped
	my ($target_ref2_mapped, $off_ref2_mapped) = qw(NA NA);
	if($design->sample_opt($sample, 'ref2_db')) {
		my $in_target = $design->get_sample_target_insert_ref2_sorted_file($sample);
		my $in_off = $design->get_sample_off_insert_ref2_sorted_file($sample);
		$target_ref2_mapped = `$samtools view $BASE_DIR/$in_target | cut -f1 | sort -u | wc -l`;
		$off_ref2_mapped = `$samtools view $BASE_DIR/$in_off | cut -f1 | sort -u | wc -l`;
		chomp $target_ref2_mapped;
		chomp $off_ref2_mapped;
	}

# get vec2 mapped
	my ($target_vec2_mapped, $off_vec2_mapped) = qw(NA NA);
	if($design->sample_opt($sample, 'vec2_db')) {
		my $in_target = $design->get_sample_target_insert_vec2_sorted_file($sample);
		my $in_off = $design->get_sample_off_insert_vec2_sorted_file($sample);
		$target_vec2_mapped = `$samtools view $BASE_DIR/$in_target | cut -f1 | sort -u | wc -l`;
		$off_vec2_mapped = `$samtools view $BASE_DIR/$in_off | cut -f1 | sort -u | wc -l`;
		chomp $target_vec2_mapped;
		chomp $off_vec2_mapped;
	}

# get target insert functional counts
	my ($target_functional_count, $target_functional_clone) = (0, 0);
	my %target_freq;
	my %target_type_freq;
	my %target_functional_type_freq;
	if($design->sample_opt($sample, 'donor_gb')) {
		my $in = $design->get_sample_target_insert_vec_summ($sample);
		open(IN, "<$BASE_DIR/$in") || die "Unable to open $BASE_DIR/$in: $!";
		<IN>; # ignore header
		while(my $line = <IN>) {
			chomp $line;
			my ($insert_id, $anno_summ, $func_clone, $type) = split(/\t/, $line);
			$target_type_freq{$type}++;
			if($func_clone > 0) {
				$target_functional_count++;
				$target_functional_clone += $func_clone;
				$target_freq{$func_clone}++;
				$target_functional_type_freq{$type}++;
			}
		}
		close(IN);
	}

# get off insert functional counts
	my ($off_functional_count, $off_functional_clone) = (0, 0);
	my %off_freq;
	my %off_type_freq;
	my %off_functional_type_freq;
	if($design->sample_opt($sample, 'donor_gb')) {
		my $in = $design->get_sample_off_insert_vec_summ($sample);
		open(IN, "<$BASE_DIR/$in") || die "Unable to open $BASE_DIR/$in: $!";
		<IN>; # ignore header
		while(my $line = <IN>) {
			chomp $line;
			my ($insert_id, $anno_summ, $func_clone, $type) = split(/\t/, $line);
			$off_type_freq{$type}++;
			if($func_clone > 0) {
				$off_functional_count++;
				$off_functional_clone += $func_clone;
				$off_freq{$func_clone}++;
				$off_functional_type_freq{$type}++;
			}
		}
		close(IN);
	}

# output
  print OUT "$sample\t$total_read\t$ref_mapped\t$ref_enrich\t$ref_target\t",
	"$target_insert\t$target_complete\t$target_incomplete\t$target_vec_mapped\t$target_nuclease_mapped\t$target_donor_mapped\t$target_trans_mapped\t$target_helper_mapped\t$target_ref2_mapped\t$target_vec2_mapped\t$target_functional_count\t$target_functional_clone\t", get_freq_str(%target_freq), "\t", get_freq_str2(%target_type_freq), "\t", get_freq_str2(%target_functional_type_freq), "\t",
	"$off_insert\t$off_complete\t$off_incomplete\t$off_vec_mapped\t$off_nuclease_mapped\t$off_donor_mapped\t$off_trans_mapped\t$off_helper_mapped\t$off_ref2_mapped\t$off_vec2_mapped\t$off_functional_count\t$off_functional_clone\t", get_freq_str(%off_freq), "\t", get_freq_str2(%off_type_freq), "\t", get_freq_str2(%off_functional_type_freq), "\n";
}

close(OUT);

sub get_freq_str {
	return 'NA' if(!@_);
	my %freq = @_;
	return join(',', map { "$_:$freq{$_}" } sort {$a <=> $b} keys %freq);
}

sub get_freq_str2 {
	return 'NA' if(!@_);
	my %freq = @_;
	return join(',', map { "$_:$freq{$_}" } sort keys %freq);
}
