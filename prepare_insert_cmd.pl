#!/bin/env perl
# Prepare bash script for extract target insert read segments and map them to vector and auxilary ref2 and vec2 genomes
our $VERSION = v1.1;
our $ENV_FILE = 'set_insert_env.sh';

use strict;
use warnings;
use lib '/project/gtplab/pipeline/Cas9_LongRead';
use Cas9OntSeqExpDesign;

my $usage = "Usage: perl $0 DESIGN-FILE BASH-OUTFILE";
my $sh_path = '/bin/bash';
my $extract_insert_script = 'extract_target_insert_seq.pl';
my $samtools = 'samtools';
my $seqret = 'seqret';

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
# prepare extract cmd
	{
		my $bed = $design->sample_opt($sample, 'target_bed');
		my $in = $design->get_sample_ref_map_target_sort_file($sample);
		my $out = $design->get_sample_target_insert_fastq($sample);
		
		my $cmd = "$SCRIPT_DIR/$extract_insert_script -t $bed -i $BASE_DIR/$in -o $BASE_DIR/$out";
		if(!(-e "$BASE_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $BASE_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare format and index cmd
	{
		my $in = $design->get_sample_target_insert_fastq($sample);
		my $out = $design->get_sample_target_insert_fasta($sample);

		my $cmd = "$seqret -sequence $BASE_DIR/$in -outseq $BASE_DIR/$out";
		$cmd .= "\n$samtools index $BASE_DIR/$out";

		if(!(-e "$BASE_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $BASE_DIR/$out already exists, won't override\n";
			$cmd =~ s/\n/\n# /sg;
			print OUT "# $cmd\n";
		}
	}

# prepare map vec cmd
	{
		my $vec_seq = $design->get_sample_vec_seq($sample);
		my $vec_map_opts = $design->sample_opt($sample, 'vec_map_opts');

		my $in = $design->get_sample_target_insert_fastq($sample);
		my $out = $design->get_sample_vec_map_file($sample);

		my $cmd = "$NGS_ALIGNER -a -x map-ont $vec_seq $BASE_DIR/$in -Y -L -t $NUM_PROC $vec_map_opts | samtools view -b -o $WORK_DIR/$out";

		if(!(-e "$WORK_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare map ref2 cmd
  if($design->sample_opt($sample, 'ref2_db')) {
		my $ref2_db = $design->sample_opt($sample, 'ref2_db');
		my $ref2_map_opts = $design->sample_opt($sample, 'ref2_map_opts');

		my $in = $design->get_sample_target_insert_fastq($sample);
		my $out = $design->get_sample_ref2_map_file($sample);

		my $cmd = "$NGS_ALIGNER -a -x map-ont $ref2_db $BASE_DIR/$in -Y -L -t $NUM_PROC $ref2_map_opts | samtools view -b -o $WORK_DIR/$out";

		if(!(-e "$WORK_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare map vec2 cmd
  if($design->sample_opt($sample, 'vec2_db')) {
		my $vec2_db = $design->sample_opt($sample, 'vec2_db');
		my $vec2_map_opts = $design->sample_opt($sample, 'vec2_map_opts');

		my $in = $design->get_sample_target_insert_fastq($sample);
		my $out = $design->get_sample_vec2_map_file($sample);

		my $cmd = "$NGS_ALIGNER -a -x map-ont $vec2_db $BASE_DIR/$in -Y -L -t $NUM_PROC $vec2_map_opts | samtools view -b -o $WORK_DIR/$out";

		if(!(-e "$WORK_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

my $min_mapQ = $design->sample_opt($sample, 'min_mapQ');
# prepare filer vec cmd
	{
		my $in = $design->get_sample_vec_map_file($sample);
		my $out = $design->get_sample_vec_filtered_file($sample);

		my $cmd = "$samtools view -q $min_mapQ -b -o $WORK_DIR/$out $WORK_DIR/$in";

		if(!(-e "$WORK_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare filer ref2 cmd
  if($design->sample_opt($sample, 'ref2_db')) {
		my $in = $design->get_sample_ref2_map_file($sample);
		my $out = $design->get_sample_ref2_filtered_file($sample);

		my $cmd = "$samtools view -q $min_mapQ -b -o $WORK_DIR/$out $WORK_DIR/$in";

		if(!(-e "$WORK_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare filer vec2 cmd
  if($design->sample_opt($sample, 'vec2_db')) {
		my $in = $design->get_sample_vec2_map_file($sample);
		my $out = $design->get_sample_vec2_filtered_file($sample);

		my $cmd = "$samtools view -q $min_mapQ -b -o $WORK_DIR/$out $WORK_DIR/$in";

		if(!(-e "$WORK_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare sort and index vec map
  {
    my $in = $design->get_sample_vec_filtered_file($sample);
    my $out = $design->get_sample_vec_sorted_file($sample);
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

# prepare sort and index ref2 map
  if($design->sample_opt($sample, 'ref2_db')) {
    my $in = $design->get_sample_ref2_filtered_file($sample);
    my $out = $design->get_sample_ref2_sorted_file($sample);
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

# prepare sort and index vec2 map
  if($design->sample_opt($sample, 'vec2_db')) {
    my $in = $design->get_sample_vec2_filtered_file($sample);
    my $out = $design->get_sample_vec2_sorted_file($sample);
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