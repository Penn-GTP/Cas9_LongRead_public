package Cas9LongReadExpDesign;
use strict;
use warnings;
use File::Basename;

our $VERSION = v1.1;
# Author: Qi Zheng
# Since: 02/01/2022

# GLOBAL options and default values
our %GLOBAL_OPTS = (
  MAX_PROC => 8,
	BASE_DIR => '.',
	WORK_DIR => 'WORK',
	SCRIPT_DIR => 'scripts',
	VEC_DIR => 'AAV_vec',
	NGS_ALIGNER => 'AAV_vec',
	FEATURE_TAG => 'label',
	MIN_COVER_RATIO => 0.9,
	ARM_KEY => '(?i:ARM)|shHDR',
 	ARM_MIN_RATIO => 0.1
);
  
# Constructor taking a filehandle or a filename
sub new {
  my $class = shift;
  my $file = shift;
  my $self = { };
  my %global_opt;
  my @sample_names;  # all samples that include these SampleIDs
  my @opt_names;  # per-sample option names
  my %sample_opt;
  # Set default global opts
	while(my ($key, $val) = each %GLOBAL_OPTS) {
		$global_opt{$key} = $val;
	}

  # read the experimental design file
  if(-f $file || -l $file) { # if $file is a file or link
		open(IN, "<$file") || die "Unable to open $file: $!";
  }
  elsif(ref $file eq 'IO') { # is a filehandle reference
	  *IN = $file; # file is a filehandle
  }
  else {
	  print STDERR "new() must take a filename or a filehandle!\n";
	  exit;
	}
  while(my $line = <IN>) {
	  chomp $line;
	  if($line =~ /^#/) {
		  if($line =~ /^## (\w+)=(\S+):/) { # global opt
			  $global_opt{$1} = $2;
		  }
		  elsif($line =~ /^## (\w+):/) { # opt name description line
			  push(@opt_names, $1);
		  }
		  else {
			  next; # ignore
		  }
	  }
	  else { # opt value line
		  my @values = split(/\t/, $line, -1);
		  if(@opt_names != @values) {
			  print STDERR "Incorrect field number at line $.: $line\n",
				  "Found ", scalar @values, " fields, but required ", scalar @opt_names, " fields\n";
			  exit;
		  }
		  my $sample = $values[0];
		  push(@sample_names, $sample);
# pair opt names with values
		  for(my $i = 0; $i < @opt_names; $i++) {
			  my $val = $values[$i];
			  $val = '' if(!defined $val);
				$sample_opt{$sample}{$opt_names[$i]} = $val;
			}
		}
	} # end each line of experiment design file

	if(-f $file) {
		close(IN);
	}

# Record variables
	%{$self->{'global_opt'}} = %global_opt;
	@{$self->{'sample_names'}} = @sample_names;
	@{$self->{'opt_names'}} = @opt_names;
	%{$self->{'sample_opt'}} = %sample_opt; 
	return bless $self, $class;
}

# Method to get global opt given opt name
sub get_global_opt {
	my ($self, $name) = @_;
	return $self->{'global_opt'}{$name};
}

# Method to get all sample_names
sub get_sample_names {
	my $self = shift;
	return @{$self->{'sample_names'}};
}

# Method to get all opt_names
sub get_opt_names {
	my $self = shift;
	return @{$self->{'opt_names'}};
}

# Method to get or set single per-sample opt value
sub sample_opt {
	my $self = shift;
	my $sample = shift;
	my $opt = shift;
	if(@_) {
		$self->{'sample_opt'}{$sample}{$opt} = shift;
	}
	return $self->{'sample_opt'}{$sample}{$opt};
}

# Method to get one or more per-sample opt values
sub get_sample_opts {
	my ($self, $sample, @opt_names) = @_;
	my @values;
	foreach my $opt (@opt_names) {
		push(@values, $self->{'sample_opt'}{$sample}{$opt});
	}
	return @values;
}

# get per-sample nuclease-vec seq file
sub get_sample_nuclease_vec_seq {
	my ($self, $sample) = @_;
	return "$sample\_nuclease_vec_seq.fasta";
}

# get per-sample nuclease-vec region file
sub get_sample_nuclease_vec_region {
	my ($self, $sample) = @_;
	return "$sample\_nuclease_vec_region.bed";
}

# get per-sample nuclease-vec annotation file
sub get_sample_nuclease_vec_anno {
	my ($self, $sample) = @_;
	return "$sample\_nuclease_vec_anno.gff3";
}

# get per-sample donor-vec seq file
sub get_sample_donor_vec_seq {
	my ($self, $sample) = @_;
	return "$sample\_donor_vec_seq.fasta";
}

# get per-sample donor-vec region file
sub get_sample_donor_vec_region {
	my ($self, $sample) = @_;
	return "$sample\_donor_vec_region.bed";
}

# get per-sample donor-vec annotation file
sub get_sample_donor_vec_anno {
	my ($self, $sample) = @_;
	return "$sample\_donor_vec_anno.gff3";
}

# get per-sample trans-vec seq file
sub get_sample_trans_vec_seq {
	my ($self, $sample) = @_;
	return "$sample\_trans_vec_seq.fasta";
}

# get per-sample trans-vec region file
sub get_sample_trans_vec_region {
	my ($self, $sample) = @_;
	return "$sample\_trans_vec_region.bed";
}

# get per-sample trans-vec annotation file
sub get_sample_trans_vec_anno {
	my ($self, $sample) = @_;
	return "$sample\_trans_vec_anno.gff3";
}

# get per-sample helper-vec seq file
sub get_sample_helper_vec_seq {
	my ($self, $sample) = @_;
	return "$sample\_helper_vec_seq.fasta";
}

# get per-sample helper-vec region file
sub get_sample_helper_vec_region {
	my ($self, $sample) = @_;
	return "$sample\_helper_vec_region.bed";
}

# get per-sample helper-vec annotation file
sub get_sample_helper_vec_anno {
	my ($self, $sample) = @_;
	return "$sample\_helper_vec_anno.gff3";
}

# get per-sample combined vec seq file
sub get_sample_vec_seq {
	my ($self, $sample) = @_;
	return "$sample\_vec_seq.fasta";
}

# get per-sample combined vec anno file
sub get_sample_vec_anno {
	my ($self, $sample) = @_;
	return "$sample\_vec_anno.gff3";
}

# get per-sample ref map file
sub get_sample_ref_map_file {
	my ($self, $sample) = @_;
	return "$sample\_ref_map.bam";
}

# get per-sample ref map filtered file
sub get_sample_ref_map_filtered_file {
	my ($self, $sample) = @_;
	return "$sample\_ref_map_filtered.bam";
}

# get per-sample ref map enrich file
sub get_sample_ref_map_enrich_file {
	my ($self, $sample) = @_;
	return "$sample\_ref_map_enrich.bam";
}

# get per-sample ref map target file
sub get_sample_ref_map_target_file {
	my ($self, $sample) = @_;
	return "$sample\_ref_map_target.bam";
}

# get per-sample ref filtered sorted file
sub get_sample_ref_map_filtered_sorted_file {
	my ($self, $sample) = @_;
	return "$sample\_ref_map_filtered_sorted.bam";
}

# get per-sample ref enrich sorted file
sub get_sample_ref_map_enrich_sorted_file {
	my ($self, $sample) = @_;
	return "$sample\_ref_map_enrich_sorted.bam";
}

# get per-sample ref target sorted file
sub get_sample_ref_map_target_sorted_file {
	my ($self, $sample) = @_;
	return "$sample\_ref_map_target_sorted.bam";
}

# get per-sample ref target insert pos
sub get_sample_target_insert_pos {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_pos.bed";
}

# get per-sample ref off insert pos
sub get_sample_off_insert_pos {
	my ($self, $sample) = @_;
	return "$sample\_off_insert_pos.bed";
}

# get per-sample ref target insert pos sorted
sub get_sample_target_insert_pos_sorted {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_pos_sorted.bed";
}

# get per-sample ref off insert pos sorted
sub get_sample_off_insert_pos_sorted {
	my ($self, $sample) = @_;
	return "$sample\_off_insert_pos_sorted.bed";
}

# get per-sample ref target insert pos merged
sub get_sample_target_insert_pos_merged {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_pos_merged.bed";
}

# get per-sample ref off insert merged
sub get_sample_off_insert_pos_merged {
	my ($self, $sample) = @_;
	return "$sample\_off_insert_pos_merged.bed";
}

# get per-sample ref off insert bed
sub get_sample_off_insert_bed {
	my ($self, $sample) = @_;
	return "$sample\_off_insert_pos.bed";
}

# get per-sample ref target insert fastq file
sub get_sample_target_insert_fastq {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_seq.fastq";
}

# get per-sample target insert fasta file
sub get_sample_target_insert_fasta {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_seq.fasta";
}

# get per-sample target insert info file
sub get_sample_target_insert_info {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_info.tsv";
}

# get per-sample ref off insert fastq file
sub get_sample_off_insert_fastq {
	my ($self, $sample) = @_;
	return "$sample\_off_insert_seq.fastq";
}

# get per-sample off insert fasta file
sub get_sample_off_insert_fasta {
	my ($self, $sample) = @_;
	return "$sample\_off_insert_seq.fasta";
}

# get per-sample off insert info file
sub get_sample_off_insert_info {
	my ($self, $sample) = @_;
	return "$sample\_off_insert_info.tsv";
}

# get per-sample target insert vec map file
sub get_sample_target_insert_vec_map_file {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_vec_map.bam";
}

# get per-sample target insert ref2 map file
sub get_sample_target_insert_ref2_map_file {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_ref2_map.bam";
}

# get per-sample target insert vec2 map file
sub get_sample_target_insert_vec2_map_file {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_vec2_map.bam";
}

# get per-sample target insert vec filtered file
sub get_sample_target_insert_vec_filtered_file {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_vec_map_filtered.bam";
}

# get per-sample target insert ref2 filtered file
sub get_sample_target_insert_ref2_filtered_file {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_ref2_map_filtered.bam";
}

# get per-sample target insert vec2 filtered file
sub get_sample_target_insert_vec2_filtered_file {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_vec2_map_filtered.bam";
}

# get per-sample target insert vec sorted file
sub get_sample_target_insert_vec_sorted_file {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_vec_map_filtered_sorted.bam";
}

# get per-sample target insert ref2 sorted file
sub get_sample_target_insert_ref2_sorted_file {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_ref2_map_filtered_sorted.bam";
}

# get per-sample target insert vec2 sorted file
sub get_sample_target_insert_vec2_sorted_file {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_vec2_map_filtered_sorted.bam";
}

# get per-sample off insert vec map file
sub get_sample_off_insert_vec_map_file {
	my ($self, $sample) = @_;
	return "$sample\_off_insert_vec_map.bam";
}

# get per-sample off insert ref2 map file
sub get_sample_off_insert_ref2_map_file {
	my ($self, $sample) = @_;
	return "$sample\_off_insert_ref2_map.bam";
}

# get per-sample off insert vec2 map file
sub get_sample_off_insert_vec2_map_file {
	my ($self, $sample) = @_;
	return "$sample\_off_insert_vec2_map.bam";
}

# get per-sample off insert vec filtered file
sub get_sample_off_insert_vec_filtered_file {
	my ($self, $sample) = @_;
	return "$sample\_off_insert_vec_map_filtered.bam";
}

# get per-sample off insert ref2 filtered file
sub get_sample_off_insert_ref2_filtered_file {
	my ($self, $sample) = @_;
	return "$sample\_off_insert_ref2_map_filtered.bam";
}

# get per-sample off insert vec2 filtered file
sub get_sample_off_insert_vec2_filtered_file {
	my ($self, $sample) = @_;
	return "$sample\_off_insert_vec2_map_filtered.bam";
}

# get per-sample off insert vec sorted file
sub get_sample_off_insert_vec_sorted_file {
	my ($self, $sample) = @_;
	return "$sample\_off_insert_vec_map_filtered_sorted.bam";
}

# get per-sample off insert ref2 sorted file
sub get_sample_off_insert_ref2_sorted_file {
	my ($self, $sample) = @_;
	return "$sample\_off_insert_ref2_map_filtered_sorted.bam";
}

# get per-sample off insert vec2 sorted file
sub get_sample_off_insert_vec2_sorted_file {
	my ($self, $sample) = @_;
	return "$sample\_off_insert_vec2_map_filtered_sorted.bam";
}

# get per-sample target insert vec inames
sub get_sample_target_insert_vec_inames {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_vec_map_inames.txt";
}

# get per-sample off insert vec inames
sub get_sample_off_insert_vec_inames {
	my ($self, $sample) = @_;
	return "$sample\_off_insert_vec_map_inames.txt";
}

# get per-sample target insert vec rnames
sub get_sample_target_insert_vec_rnames {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_vec_map_rnames.txt";
}

# get per-sample off insert vec rnames
sub get_sample_off_insert_vec_rnames {
	my ($self, $sample) = @_;
	return "$sample\_off_insert_vec_map_rnames.txt";
}

# get per-sample target insert vec fasta
sub get_sample_target_insert_vec_fasta {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_vec_map_seq.fasta";
}

# get per-sample off insert vec fasta
sub get_sample_off_insert_vec_fasta {
	my ($self, $sample) = @_;
	return "$sample\_off_insert_vec_map_seq.fasta";
}

# get per-sample target insert ref sorted file
sub get_sample_target_insert_ref_sorted_file {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_ref_map_sorted.bam";
}

# get per-sample off insert ref sorted file
sub get_sample_off_insert_ref_sorted_file {
	my ($self, $sample) = @_;
	return "$sample\_off_insert_ref_map_sorted.bam";
}

# get per-sample target_genomic ref sorted file
sub get_sample_target_genomic_ref_sorted_file {
	my ($self, $sample) = @_;
	return "$sample\_target_genomic_ref_map_sorted.bam";
}

# get per-sample target_genomic ref map seq
sub get_sample_target_genomic_ref_map_seq {
	my ($self, $sample) = @_;
	return "$sample\_target_genomic_ref_map_seq.fasta";
}

# get per-sample target insert size distribution figure
sub get_sample_target_insert_size_distrib {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_size_distrib.pdf";
}

# get per-sample off insert size distribution figure
sub get_sample_off_insert_size_distrib {
	my ($self, $sample) = @_;
	return "$sample\_off_insert_size_distrib.pdf";
}

# get per-sample target insert vec anno
sub get_sample_target_insert_vec_anno {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_vec_anno.gff3";
}

# get per-sample target insert nuclease_vec anno
sub get_sample_target_insert_nuclease_vec_anno {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_nuclease_vec_anno.gff3";
}

# get per-sample target insert donor_vec anno
sub get_sample_target_insert_donor_vec_anno {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_donor_vec_anno.gff3";
}

# get per-sample target insert trans_vec anno
sub get_sample_target_insert_trans_vec_anno {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_trans_vec_anno.gff3";
}

# get per-sample target insert helper_vec anno
sub get_sample_target_insert_helper_vec_anno {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_helper_vec_anno.gff3";
}

# get per-sample target insert ref2 anno file
sub get_sample_target_insert_ref2_anno {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_ref2_anno.gff3";
}

# get per-sample target insert vec2 anno file
sub get_sample_target_insert_vec2_anno {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_vec2_anno.gff3";
}

# get per-sample off insert vec anno
sub get_sample_off_insert_vec_anno {
	my ($self, $sample) = @_;
	return "$sample\_off_insert_vec_anno.gff3";
}

# get per-sample off insert nuclease_vec anno
sub get_sample_off_insert_nuclease_vec_anno {
	my ($self, $sample) = @_;
	return "$sample\_off_insert_nuclease_vec_anno.gff3";
}

# get per-sample off insert donor_vec anno
sub get_sample_off_insert_donor_vec_anno {
	my ($self, $sample) = @_;
	return "$sample\_off_insert_donor_vec_anno.gff3";
}

# get per-sample off insert trans_vec anno
sub get_sample_off_insert_trans_vec_anno {
	my ($self, $sample) = @_;
	return "$sample\_off_insert_trans_vec_anno.gff3";
}

# get per-sample off insert helper_vec anno
sub get_sample_off_insert_helper_vec_anno {
	my ($self, $sample) = @_;
	return "$sample\_off_insert_helper_vec_anno.gff3";
}

# get per-sample off insert ref2 anno file
sub get_sample_off_insert_ref2_anno {
	my ($self, $sample) = @_;
	return "$sample\_off_insert_ref2_anno.gff3";
}

# get per-sample off insert vec2 anno file
sub get_sample_off_insert_vec2_anno {
	my ($self, $sample) = @_;
	return "$sample\_off_insert_vec2_anno.gff3";
}

# get per-sample target genomic ref anno
sub get_sample_target_genomic_ref_anno {
	my ($self, $sample) = @_;
	return "$sample\_target_genomic_ref_anno.gff3";
}

# get per-sample insert nuclease_vec feature count
sub get_sample_insert_nuclease_vec_count {
	my ($self, $sample) = @_;
	return "$sample\_insert_nuclease_vec_feature_count.tsv";
}

# get per-sample insert donor_vec feature count
sub get_sample_insert_donor_vec_count {
	my ($self, $sample) = @_;
	return "$sample\_insert_donor_vec_feature_count.tsv";
}

# get per-sample target insert vec summ
sub get_sample_target_insert_vec_summ {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_vec_summ.tsv";
}

# get per-sample off insert vec summ
sub get_sample_off_insert_vec_summ {
	my ($self, $sample) = @_;
	return "$sample\_off_insert_vec_summ.tsv";
}

# get per-exp stats
sub get_exp_stats_file {
	my ($self, $exp_file) = @_;
	my $stats_file = basename($exp_file, qw(.conf .txt .tsv));
	return "$stats_file\_sample_stats.tsv";
}

1;
