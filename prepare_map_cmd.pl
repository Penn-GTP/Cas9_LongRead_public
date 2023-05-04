#!/bin/env perl
# Prepare bash script for mapping Cas9 LongReads to ref genome and extract inserts
our $VERSION = 'v1.1.1';
our $ENV_FILE = 'set_map_env.sh';

use strict;
use warnings;
use lib '/project/gtplab/pipeline/Cas9_LongRead';
use Cas9LongReadExpDesign;

my $usage = "Usage: perl $0 DESIGN-FILE BASH-OUTFILE";
my $sh_path = '/bin/bash';
my $samtools = 'samtools';
my $cmd = "$0 " . join(" ", @ARGV);

my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $design = new Cas9LongReadExpDesign($infile);
my $NUM_PROC = $design->get_global_opt('NUM_PROC');
my $BASE_DIR = $design->get_global_opt('BASE_DIR');
my $SCRIPT_DIR = $design->get_global_opt('SCRIPT_DIR');
my $VEC_DIR = $design->get_global_opt('VEC_DIR');
my $WORK_DIR = $design->get_global_opt('WORK_DIR');
my $NGS_ALIGNER = $design->get_global_opt('NGS_ALIGNER');

my $DEFAULT_TECH = 'ont';

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
print OUT qq(# CMD:"$cmd"\n# VER:$VERSION\n);
# set env
print OUT "source $SCRIPT_DIR/$ENV_FILE\n\n";

foreach my $sample ($design->get_sample_names()) {
# prepare map cmd
	{
		my $ref_db = $design->sample_opt($sample, 'ref_db');
		my $tech = $design->sample_opt($sample, 'longread_tech') ? $design->sample_opt($sample, 'longread_tech') : $DEFAULT_TECH;
		my $ref_map_opts = $design->sample_opt($sample, 'ref_map_opts');
		my $in = $design->sample_opt($sample, 'read_fastq');
		my $out = $design->get_sample_ref_map_file($sample);

		my $cmd = "$NGS_ALIGNER -a -x map-$tech $ref_db $in -Y -L -t $NUM_PROC $ref_map_opts | samtools view -b -o $WORK_DIR/$out";

		if(!(-e "$WORK_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare filtered alignment
  {
    my $in = $design->get_sample_ref_map_file($sample);
    my $out = $design->get_sample_ref_map_filtered_file($sample);
		my $min_mapQ = $design->sample_opt($sample, 'min_mapQ');

    my $cmd = "$samtools view -q $min_mapQ -F 0x4 -F 0x100 -b -o $WORK_DIR/$out $WORK_DIR/$in"; # no secondary alignments

    if(!-e "$WORK_DIR/$out") {
      print OUT "$cmd\n";
    }
    else {
			print STDERR "Warning: $WORK_DIR/$out already exists, won't override\n";
      print OUT "# $cmd\n";
    }
  }

# prepare extract enrich alignment
  {
    my $in = $design->get_sample_ref_map_filtered_file($sample);
		my $bed = $design->sample_opt($sample, 'enrich_bed');
    my $out = $design->get_sample_ref_map_enrich_file($sample);
		my $opts = $bed && -e $bed ? "-L $bed" : ""; # test whether enrich_bed exists

    my $cmd = "$samtools view $opts -b -o $WORK_DIR/$out $WORK_DIR/$in";

    if(!-e "$WORK_DIR/$out") {
      print OUT "$cmd\n";
    }
    else {
			print STDERR "Warning: $WORK_DIR/$out already exists, won't override\n";
      print OUT "# $cmd\n";
    }
  }

# prepare extract target alignment
  {
    my $in = $design->get_sample_ref_map_enrich_file($sample);
		my $bed = $design->sample_opt($sample, 'target_bed');
    my $out = $design->get_sample_ref_map_target_file($sample);
    my $cmd = "$samtools view -L $bed -b -o $WORK_DIR/$out $WORK_DIR/$in";

    if(!-e "$WORK_DIR/$out") {
      print OUT "$cmd\n";
    }
    else {
			print STDERR "Warning: $WORK_DIR/$out already exists, won't override\n";
      print OUT "# $cmd\n";
    }
  }

# sorted and index filtered file
  {
    my $in = $design->get_sample_ref_map_filtered_file($sample);
    my $out = $design->get_sample_ref_map_filtered_sorted_file($sample);
    my $cmd = "$samtools sort $WORK_DIR/$in -o $BASE_DIR/$out";
		$cmd .= "\n$samtools index $BASE_DIR/$out";

    if(!-e "$BASE_DIR/$out") {
      print OUT "$cmd\n";
    }
    else {
			print STDERR "Warning: $BASE_DIR/$out already exists, won't override\n";
			$cmd =~ s/\n/\n# /sg;
      print OUT "# $cmd\n";
    }
  }

# sorted and index enrich file
  {
    my $in = $design->get_sample_ref_map_enrich_file($sample);
    my $out = $design->get_sample_ref_map_enrich_sorted_file($sample);
    my $cmd = "$samtools sort $WORK_DIR/$in -o $BASE_DIR/$out";
		$cmd .= "\n$samtools index $BASE_DIR/$out";

    if(!-e "$BASE_DIR/$out") {
      print OUT "$cmd\n";
    }
    else {
			print STDERR "Warning: $BASE_DIR/$out already exists, won't override\n";
			$cmd =~ s/\n/\n# /sg;
      print OUT "# $cmd\n";
    }
  }

# sorted and index target file
  {
    my $in = $design->get_sample_ref_map_target_file($sample);
    my $out = $design->get_sample_ref_map_target_sorted_file($sample);
    my $cmd = "$samtools sort $WORK_DIR/$in -o $BASE_DIR/$out";
		$cmd .= "\n$samtools index $BASE_DIR/$out";

    if(!-e "$BASE_DIR/$out") {
      print OUT "$cmd\n";
    }
    else {
			print STDERR "Warning: $BASE_DIR/$out already exists, won't override\n";
			$cmd =~ s/\n/\n# /sg;
      print OUT "# $cmd\n";
    }
  }

  print OUT "\n";
}

close(OUT);
# change to exacutable
chmod 0750, $outfile;
