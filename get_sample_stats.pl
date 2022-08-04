#!/bin/env perl
# Prepare bash script for generating per-sample processing stats
our $VERSION = v1.1;

use strict;
use warnings;
use lib '/project/gtplab/pipeline/Cas9_LongRead';
use Cas9LongReadExpDesign;

my $usage = "Usage: perl $0 DESIGN-FILE OUTFILE";
#my $sh_path = '/bin/bash';
my $samtools = 'samtools';
my @headers = qw(sample_name total_read ref_mapped ref_enrich ref_target
target_insert target_insert_complete target_insert_incomplete target_insert_vec_mapped target_insert_nuclease_mapped target_insert_donor_mapped target_insert_trans_mapped target_insert_helper_mapped target_insert_ref2_mapped target_insert_vec2_mapped target_insert_functional_donor_basic target_insert_functional_donor_full
off_insert off_insert_complete off_insert_incomplete off_insert_vec_mapped off_insert_nuclease_mapped off_insert_donor_mapped off_insert_trans_mapped off_insert_helper_mapped off_insert_ref2_mapped off_insert_vec2_mapped off_insert_functional_donor_basic off_insert_functional_donor_full);

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
	if($design->get_sample_opts($sample, 'ref2_db')) {
		my $in_target = $design->get_sample_target_insert_ref2_sorted_file($sample);
		my $in_off = $design->get_sample_off_insert_ref2_sorted_file($sample);
		$target_ref2_mapped = `$samtools view $BASE_DIR/$in_target | cut -f1 | sort -u | wc -l`;
		$off_ref2_mapped = `$samtools view $BASE_DIR/$in_off | cut -f1 | sort -u | wc -l`;
		chomp $target_ref2_mapped;
		chomp $off_ref2_mapped;
	}

# get vec2 mapped
	my ($target_vec2_mapped, $off_vec2_mapped) = qw(NA NA);
	if($design->get_sample_opts($sample, 'vec2_db')) {
		my $in_target = $design->get_sample_target_insert_vec2_sorted_file($sample);
		my $in_off = $design->get_sample_off_insert_vec2_sorted_file($sample);
		$target_vec2_mapped = `$samtools view $BASE_DIR/$in_target | cut -f1 | sort -u | wc -l`;
		$off_vec2_mapped = `$samtools view $BASE_DIR/$in_off | cut -f1 | sort -u | wc -l`;
		chomp $target_vec2_mapped;
		chomp $off_vec2_mapped;
	}

# get functional donor counts
	my ($target_functional_donor_basic, $off_functional_donor_basic, $target_functional_donor_full, $off_functional_donor_full) = (0, 0, 0, 0);
	if($design->sample_opt($sample, 'donor_gb')) {
		my @feat_basic = split(/\|/, $design->sample_opt($sample, 'functional_donor_feature_basic'));
		my @feat_full = split(/\|/, $design->sample_opt($sample, 'functional_donor_feature_full'));
    {
		  open(GFF, "<", $design->get_sample_target_insert_donor_vec_anno($sample));
			my %name2feat = get_anno_summ(\*GFF);
			foreach my $feats (values %name2feat) {
				my @feats_rev = reverse(@$feats);
				if(is_sub_array($feats, \@feat_basic) || is_sub_array(\@feats_rev, \@feat_basic)) {
					$target_functional_donor_basic++;
				}
				if(is_sub_array($feats, \@feat_full) || is_sub_array(\@feats_rev, \@feat_full)) {
					$target_functional_donor_full++;
				}
			}
			close(GFF);
    }

    {
		  open(GFF, "<", $design->get_sample_off_insert_donor_vec_anno($sample));
			my %name2feat = get_anno_summ(\*GFF);
			foreach my $feats (values %name2feat) {
				my @feats_rev = reverse(@$feats);
				if(is_sub_array($feats, \@feat_basic) || is_sub_array(\@feats_rev, \@feat_basic)) {
					$off_functional_donor_basic++;
				}
				if(is_sub_array($feats, \@feat_full) || is_sub_array(\@feats_rev, \@feat_full)) {
					$off_functional_donor_full++;
				}
			}
			close(GFF);
    }
	}

# output
  print OUT "$sample\t$total_read\t$ref_mapped\t$ref_enrich\t$ref_target\t",
	"$target_insert\t$target_complete\t$target_incomplete\t$target_vec_mapped\t$target_nuclease_mapped\t$target_donor_mapped\t$target_trans_mapped\t$target_helper_mapped\t$target_ref2_mapped\t$target_vec2_mapped\t$target_functional_donor_basic\t$target_functional_donor_full\t",
	"$off_insert\t$off_complete\t$off_incomplete\t$off_vec_mapped\t$off_nuclease_mapped\t$off_donor_mapped\t$off_trans_mapped\t$off_helper_mapped\t$off_ref2_mapped\t$off_vec2_mapped\t$off_functional_donor_basic\t$off_functional_donor_full\n";
}

close(OUT);

# subroutine definitions
sub get_anno_summ {
	my $fh = shift;
	my %name2feat;
	while(my $line = <$fh>) {
		chomp $line;
		next if($line =~ /^#/);
		my ($name, $src, $type, $start, $end, $score, $strand, $frame, $attrs) = split(/\t/, $line);
    foreach my $attr (split(/;/, $attrs)) {
			my ($tag, $val) = split(/=/, $attr);
			if($tag eq $FEATURE_TAG) {
				push(@{$name2feat{$name}}, $val);
			}
		}
	}
	return %name2feat;
}

sub is_sub_array {
	my ($A, $B) = @_;
	my $n = @$A;
	my $m = @$B;
	my ($i, $j) = (0, 0);
	while($i < $n && $j < $m) {
		if($A->[$i] eq $B->[$j]) {
			$i++;
			$j++;

			if($j == $m) {
				return 1;
			}
		}
		else {
			$i = $i - $j + 1;
			$j = 0;
		}
	}

	return 0;
}

