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
	NGS_ALIGNER => 'AAV_vec'
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

# get per-sample ref enrich sort file
sub get_sample_ref_map_enrich_sort_file {
	my ($self, $sample) = @_;
	return "$sample\_ref_map_enrich_sorted.bam";
}

# get per-sample ref target sort file
sub get_sample_ref_map_target_sort_file {
	my ($self, $sample) = @_;
	return "$sample\_ref_map_target_sorted.bam";
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

# get per-sample vec map file
sub get_sample_vec_map_file {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_vec_map.bam";
}

# get per-sample ref2 map file
sub get_sample_ref2_map_file {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_ref2_map.bam";
}

# get per-sample vec2 map file
sub get_sample_vec2_map_file {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_vec2_map.bam";
}

# get per-sample vec filtered file
sub get_sample_vec_filtered_file {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_vec_map_filtered.bam";
}

# get per-sample ref2 filtered file
sub get_sample_ref2_filtered_file {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_ref2_map_filtered.bam";
}

# get per-sample vec2 filtered file
sub get_sample_vec2_filtered_file {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_vec2_map_filtered.bam";
}

# get per-sample vec sorted file
sub get_sample_vec_sorted_file {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_vec_map_filtered_sorted.bam";
}

# get per-sample ref2 sorted file
sub get_sample_ref2_sorted_file {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_ref2_map_filtered_sorted.bam";
}

# get per-sample vec2 sorted file
sub get_sample_vec2_sorted_file {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_vec2_map_filtered_sorted.bam";
}

# get per-sample target insert size distribution figure
sub get_sample_target_insert_size_distrib {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_size_distrib.pdf";
}

# get per-sample vec sorted bed
sub get_sample_vec_sorted_bed {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_vec_map_filtered_sorted.bed";
}

# get per-sample ref2 sorted bed
sub get_sample_ref2_sorted_bed {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_ref2_map_filtered_sorted.bed";
}

# get per-sample vec2 sorted bed
sub get_sample_vec2_sorted_bed {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_vec2_map_filtered_sorted.bed";
}

# get per-sample insert nuclease_vec anno
sub get_sample_insert_nuclease_vec_anno {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_nuclease_vec_anno.gff3";
}

# get per-sample donor_vec anno
sub get_sample_insert_donor_vec_anno {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_donor_vec_anno.gff3";
}

# get per-sample insert trans_vec anno
sub get_sample_insert_trans_vec_anno {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_trans_vec_anno.gff3";
}

# get per-sample helper_vec anno
sub get_sample_insert_helper_vec_anno {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_helper_vec_anno.gff3";
}

# get per-sample insert ref2 anno file
sub get_sample_insert_ref2_anno {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_ref2_anno.gff3";
}

# get per-sample insert vec2 anno file
sub get_sample_insert_vec2_anno {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_vec2_anno.gff3";
}

# get per-sample insert nuclease_vec feature count
sub get_sample_insert_nuclease_vec_count {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_nuclease_vec_feature_count.tsv";
}

# get per-sample insert donor_vec feature count
sub get_sample_insert_donor_vec_count {
	my ($self, $sample) = @_;
	return "$sample\_target_insert_donor_vec_feature_count.tsv";
}

# get per-sample stats
sub get_sample_stats_file {
	my ($self, $sample) = @_;
	return "$sample\_sample_stats.tsv";
}

1;