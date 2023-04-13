#!/bin/env perl
# Prepare bash script for extract target insert read segments and map them to vector and auxilary ref2 and vec2 genomes
our $VERSION = v1.1;
our $ENV_FILE = 'set_insert_env.sh';

use strict;
use warnings;
use lib '/project/gtplab/pipeline/Cas9_LongRead';
use Cas9LongReadExpDesign;

my $usage = "Usage: perl $0 DESIGN-FILE BASH-OUTFILE";
my $sh_path = '/bin/bash';
my $target_pos_script = 'get_target_insert_pos.pl';
my $off_pos_script = 'get_off_insert_pos.pl';
my $extract_insert_script = 'extract_insert.pl';
my $filter_fasta_script = 'filter_fasta_file.pl';
my $samtools = 'samtools';
my $bedtools = 'bedtools';
my $picard = 'picard.jar';

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
my $DEFAULT_MIN_INSERT = 20;
my $DEFAULT_MAX_DIST = 520;

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
	my $tech = $design->sample_opt($sample, 'longread_tech') ? $design->sample_opt($sample, 'longread_tech') : $DEFAULT_TECH;
	my $min_insert = $design->sample_opt($sample, 'min_insert') ? $design->sample_opt($sample, 'min_insert') : $DEFAULT_MIN_INSERT;
	my $max_dist = $design->sample_opt($sample, 'max_dist') ? $design->sample_opt($sample, 'max_dist') : $DEFAULT_MAX_DIST;
# prepare get target insert pos cmd
	{
		my $bed = $design->sample_opt($sample, 'target_bed');
		my $in = $design->get_sample_ref_map_target_sorted_file($sample);
		my $out = $design->get_sample_target_insert_pos($sample);
		my $sorted = $design->get_sample_target_insert_pos_sorted($sample);
		my $merged = $design->get_sample_target_insert_pos_merged($sample);
		
		my $cmd = "$SCRIPT_DIR/$target_pos_script -t $bed -i $BASE_DIR/$in -o $WORK_DIR/$out --min-insert $min_insert --max-dist $max_dist";
		$cmd .= "\nif [ -s $WORK_DIR/$out ];";
		$cmd .= "\nthen $bedtools sort -i $WORK_DIR/$out > $WORK_DIR/$sorted;";
		$cmd .= "\n$bedtools merge -s -d $max_dist -i $WORK_DIR/$sorted -c 4,5,6 -o collapse,sum,distinct > $WORK_DIR/$merged;";
		$cmd .= "\nelse cp $WORK_DIR/$out $WORK_DIR/$sorted;";
		$cmd .= "\ncp $WORK_DIR/$sorted $WORK_DIR/$merged;";
		$cmd .= "\nfi";

		if(!(-e "$WORK_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out already exists, won't override\n";
			$cmd =~ s/\n/\n# /sg;
			print OUT "# $cmd\n";
		}
	}

# prepare get off insert pos cmd
	{
		my $bed = $design->sample_opt($sample, 'target_bed');
		my $in = $design->get_sample_ref_map_filtered_sorted_file($sample);
		my $out = $design->get_sample_off_insert_pos($sample);
		my $sorted = $design->get_sample_off_insert_pos_sorted($sample);
		my $merged = $design->get_sample_off_insert_pos_merged($sample);
		
		my $cmd = "$SCRIPT_DIR/$off_pos_script -t $bed -i $BASE_DIR/$in -o $WORK_DIR/$out --min-insert $min_insert --max-dist $max_dist";
		$cmd .= "\nif [ -s $WORK_DIR/$out ];";
		$cmd .= "\nthen $bedtools sort -i $WORK_DIR/$out > $WORK_DIR/$sorted;";
		$cmd .= "\n$bedtools merge -s -d $max_dist -i $WORK_DIR/$sorted -c 4,5,6 -o collapse,sum,distinct > $WORK_DIR/$merged;";
		$cmd .= "\nelse cp $WORK_DIR/$out $WORK_DIR/$sorted;";
		$cmd .= "\ncp $WORK_DIR/$sorted $WORK_DIR/$merged;";
		$cmd .= "\nfi";

		if(!(-e "$WORK_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out already exists, won't override\n";
			$cmd =~ s/\n/\n# /sg;
			print OUT "# $cmd\n";
		}
	}

# prepare extract target and index cmd
	{
		my $in = $design->get_sample_target_insert_pos_merged($sample);
		my $read = $design->sample_opt($sample, 'read_fastq');
		my $fq_out = $design->get_sample_target_insert_fastq($sample);
		my $fa_out = $design->get_sample_target_insert_fasta($sample);
		my $info_out = $design->get_sample_target_insert_info($sample);
		
		my $cmd = "$SCRIPT_DIR/$extract_insert_script -i $WORK_DIR/$in -s $read -fq $BASE_DIR/$fq_out -fa $BASE_DIR/$fa_out -info $BASE_DIR/$info_out";
		$cmd .= "\n$samtools faidx $BASE_DIR/$fa_out";

		if(!(-e "$BASE_DIR/$fq_out" && -e "$BASE_DIR/$fa_out" && -e "$BASE_DIR/$info_out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $BASE_DIR/$fq_out, $BASE_DIR/$fa_out and $BASE_DIR/$info_out already exists, won't override\n";
			$cmd =~ s/\n/\n# /sg;
			print OUT "# $cmd\n";
		}
	}

# prepare extract off and index cmd
	{
		my $in = $design->get_sample_off_insert_pos_merged($sample);
		my $read = $design->sample_opt($sample, 'read_fastq');
		my $fq_out = $design->get_sample_off_insert_fastq($sample);
		my $fa_out = $design->get_sample_off_insert_fasta($sample);
		my $info_out = $design->get_sample_off_insert_info($sample);
		
		my $cmd = "$SCRIPT_DIR/$extract_insert_script -i $WORK_DIR/$in -s $read -fq $BASE_DIR/$fq_out -fa $BASE_DIR/$fa_out -info $BASE_DIR/$info_out";
		$cmd .= "\n$samtools faidx $BASE_DIR/$fa_out";

		if(!(-e "$BASE_DIR/$fq_out" && -e "$BASE_DIR/$fa_out" && -e "$BASE_DIR/$info_out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $BASE_DIR/$fq_out, $BASE_DIR/$fa_out and $BASE_DIR/$info_out already exists, won't override\n";
			$cmd =~ s/\n/\n# /sg;
			print OUT "# $cmd\n";
		}
	}

# prepare map target insert vec cmd
	{
		my $vec_db = $design->get_sample_vec_seq($sample);
		my $vec_map_opts = $design->sample_opt($sample, 'vec_map_opts');

		my $in = $design->get_sample_target_insert_fastq($sample);
		my $out = $design->get_sample_target_insert_vec_map_file($sample);

		my $cmd = "$NGS_ALIGNER -a -x map-$tech $VEC_DIR/$vec_db $BASE_DIR/$in -Y -L -t $NUM_PROC $vec_map_opts | samtools view -b -o $WORK_DIR/$out";

		if(!(-e "$WORK_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare map target insert ref2 cmd
  if($design->sample_opt($sample, 'ref2_db')) {
		my $ref2_db = $design->sample_opt($sample, 'ref2_db');
		my $ref2_map_opts = $design->sample_opt($sample, 'ref2_map_opts');

		my $in = $design->get_sample_target_insert_fastq($sample);
		my $out = $design->get_sample_target_insert_ref2_map_file($sample);

		my $cmd = "$NGS_ALIGNER -a -x map-$tech $ref2_db $BASE_DIR/$in -Y -L -t $NUM_PROC $ref2_map_opts | samtools view -b -o $WORK_DIR/$out";

		if(!(-e "$WORK_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare map target insert vec2 cmd
  if($design->sample_opt($sample, 'vec2_db')) {
		my $vec2_db = $design->sample_opt($sample, 'vec2_db');
		my $vec2_map_opts = $design->sample_opt($sample, 'vec2_map_opts');

		my $in = $design->get_sample_target_insert_fastq($sample);
		my $out = $design->get_sample_target_insert_vec2_map_file($sample);

		my $cmd = "$NGS_ALIGNER -a -x map-$tech $vec2_db $BASE_DIR/$in -Y -L -t $NUM_PROC $vec2_map_opts | samtools view -b -o $WORK_DIR/$out";

		if(!(-e "$WORK_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare map off insert vec cmd
	{
		my $vec_db = $design->get_sample_vec_seq($sample);
		my $vec_map_opts = $design->sample_opt($sample, 'vec_map_opts');

		my $in = $design->get_sample_off_insert_fastq($sample);
		my $out = $design->get_sample_off_insert_vec_map_file($sample);

		my $cmd = "$NGS_ALIGNER -a -x map-$tech $VEC_DIR/$vec_db $BASE_DIR/$in -Y -L -t $NUM_PROC $vec_map_opts | samtools view -b -o $WORK_DIR/$out";

		if(!(-e "$WORK_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare map off insert ref2 cmd
  if($design->sample_opt($sample, 'ref2_db')) {
		my $ref2_db = $design->sample_opt($sample, 'ref2_db');
		my $ref2_map_opts = $design->sample_opt($sample, 'ref2_map_opts');

		my $in = $design->get_sample_off_insert_fastq($sample);
		my $out = $design->get_sample_off_insert_ref2_map_file($sample);

		my $cmd = "$NGS_ALIGNER -a -x map-$tech $ref2_db $BASE_DIR/$in -Y -L -t $NUM_PROC $ref2_map_opts | samtools view -b -o $WORK_DIR/$out";

		if(!(-e "$WORK_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare map off insert vec2 cmd
  if($design->sample_opt($sample, 'vec2_db')) {
		my $vec2_db = $design->sample_opt($sample, 'vec2_db');
		my $vec2_map_opts = $design->sample_opt($sample, 'vec2_map_opts');

		my $in = $design->get_sample_off_insert_fastq($sample);
		my $out = $design->get_sample_off_insert_vec2_map_file($sample);

		my $cmd = "$NGS_ALIGNER -a -x map-$tech $vec2_db $BASE_DIR/$in -Y -L -t $NUM_PROC $vec2_map_opts | samtools view -b -o $WORK_DIR/$out";

		if(!(-e "$WORK_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

my $min_mapQ = $design->sample_opt($sample, 'min_mapQ');

# prepare target insert filtered vec cmd
	{
		my $in = $design->get_sample_target_insert_vec_map_file($sample);
		my $out = $design->get_sample_target_insert_vec_filtered_file($sample);

		my $cmd = "$samtools view -q $min_mapQ -b -o $WORK_DIR/$out $WORK_DIR/$in";

		if(!(-e "$WORK_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare target insert filtered ref2 cmd
  if($design->sample_opt($sample, 'ref2_db')) {
		my $in = $design->get_sample_target_insert_ref2_map_file($sample);
		my $out = $design->get_sample_target_insert_ref2_filtered_file($sample);

		my $cmd = "$samtools view -q $min_mapQ -b -o $WORK_DIR/$out $WORK_DIR/$in";

		if(!(-e "$WORK_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare target insert filtered vec2 cmd
  if($design->sample_opt($sample, 'vec2_db')) {
		my $in = $design->get_sample_target_insert_vec2_map_file($sample);
		my $out = $design->get_sample_target_insert_vec2_filtered_file($sample);

		my $cmd = "$samtools view -q $min_mapQ -b -o $WORK_DIR/$out $WORK_DIR/$in";

		if(!(-e "$WORK_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare target insert sorted and index vec map
  {
    my $in = $design->get_sample_target_insert_vec_filtered_file($sample);
    my $out = $design->get_sample_target_insert_vec_sorted_file($sample);
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

# prepare target insert sorted and index ref2 map
  if($design->sample_opt($sample, 'ref2_db')) {
    my $in = $design->get_sample_target_insert_ref2_filtered_file($sample);
    my $out = $design->get_sample_target_insert_ref2_sorted_file($sample);
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

# prepare target insert sorted and index vec2 map
  if($design->sample_opt($sample, 'vec2_db')) {
    my $in = $design->get_sample_target_insert_vec2_filtered_file($sample);
    my $out = $design->get_sample_target_insert_vec2_sorted_file($sample);
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

# prepare off insert filtered vec cmd
	{
		my $in = $design->get_sample_off_insert_vec_map_file($sample);
		my $out = $design->get_sample_off_insert_vec_filtered_file($sample);

		my $cmd = "$samtools view -q $min_mapQ -b -o $WORK_DIR/$out $WORK_DIR/$in";

		if(!(-e "$WORK_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare off insert filtered ref2 cmd
  if($design->sample_opt($sample, 'ref2_db')) {
		my $in = $design->get_sample_off_insert_ref2_map_file($sample);
		my $out = $design->get_sample_off_insert_ref2_filtered_file($sample);

		my $cmd = "$samtools view -q $min_mapQ -b -o $WORK_DIR/$out $WORK_DIR/$in";

		if(!(-e "$WORK_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare off insert filtered vec2 cmd
  if($design->sample_opt($sample, 'vec2_db')) {
		my $in = $design->get_sample_off_insert_vec2_map_file($sample);
		my $out = $design->get_sample_off_insert_vec2_filtered_file($sample);

		my $cmd = "$samtools view -q $min_mapQ -b -o $WORK_DIR/$out $WORK_DIR/$in";

		if(!(-e "$WORK_DIR/$out")) {
			print OUT "$cmd\n";
		}
		else {
			print STDERR "Warning: $WORK_DIR/$out already exists, won't override\n";
			print OUT "# $cmd\n";
		}
	}

# prepare off insert sorted and index vec map
  {
    my $in = $design->get_sample_off_insert_vec_filtered_file($sample);
    my $out = $design->get_sample_off_insert_vec_sorted_file($sample);
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

# prepare off insert sorted and index ref2 map
  if($design->sample_opt($sample, 'ref2_db')) {
    my $in = $design->get_sample_off_insert_ref2_filtered_file($sample);
    my $out = $design->get_sample_off_insert_ref2_sorted_file($sample);
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

# prepare off insert sorted and index vec2 map
  if($design->sample_opt($sample, 'vec2_db')) {
    my $in = $design->get_sample_off_insert_vec2_filtered_file($sample);
    my $out = $design->get_sample_off_insert_vec2_sorted_file($sample);
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

# prepare target insert vec inames
  {
    my $in = $design->get_sample_target_insert_vec_sorted_file($sample);
    my $out = $design->get_sample_target_insert_vec_inames($sample);
    my $cmd = "$samtools view $BASE_DIR/$in | cut -f1 | sort -u > $WORK_DIR/$out";

    if(!-e "$WORK_DIR/$out") {
      print OUT "$cmd\n";
    }
    else {
			print STDERR "Warning: $WORK_DIR/$out already exists, won't override\n";
      print OUT "# $cmd\n";
    }
  }

# prepare off insert vec inames
  {
    my $in = $design->get_sample_off_insert_vec_sorted_file($sample);
    my $out = $design->get_sample_off_insert_vec_inames($sample);
    my $cmd = "$samtools view $BASE_DIR/$in | cut -f1 | sort -u > $WORK_DIR/$out";

    if(!-e "$WORK_DIR/$out") {
      print OUT "$cmd\n";
    }
    else {
			print STDERR "Warning: $WORK_DIR/$out already exists, won't override\n";
      print OUT "# $cmd\n";
    }
  }

# prepare target insert vec rnames
  {
    my $in = $design->get_sample_target_insert_vec_sorted_file($sample);
    my $out = $design->get_sample_target_insert_vec_rnames($sample);
    my $cmd = "$samtools view $BASE_DIR/$in | cut -f1 | perl -p -e 's/:.*//' | sort -u > $WORK_DIR/$out";

    if(!-e "$WORK_DIR/$out") {
      print OUT "$cmd\n";
    }
    else {
			print STDERR "Warning: $WORK_DIR/$out already exists, won't override\n";
      print OUT "# $cmd\n";
    }
  }

# prepare off insert vec rnames
  {
    my $in = $design->get_sample_off_insert_vec_sorted_file($sample);
    my $out = $design->get_sample_off_insert_vec_rnames($sample);
    my $cmd = "$samtools view $BASE_DIR/$in | cut -f1 | perl -p -e 's/:.*//' | sort -u > $WORK_DIR/$out";

    if(!-e "$WORK_DIR/$out") {
      print OUT "$cmd\n";
    }
    else {
			print STDERR "Warning: $WORK_DIR/$out already exists, won't override\n";
      print OUT "# $cmd\n";
    }
  }

# prepare target insert vec fasta and index
  {
    my $in = $design->get_sample_target_insert_fasta($sample);
    my $list = $design->get_sample_target_insert_vec_inames($sample);
    my $out = $design->get_sample_target_insert_vec_fasta($sample);
    my $cmd = "$filter_fasta_script -i $BASE_DIR/$in -l $WORK_DIR/$list -o $BASE_DIR/$out";
		$cmd .= "\n$samtools faidx $BASE_DIR/$out";

    if(!-e "$BASE_DIR/$out") {
      print OUT "$cmd\n";
    }
    else {
			print STDERR "Warning: $BASE_DIR/$out already exists, won't override\n";
			$cmd =~ s/\n/\n# /sg;
      print OUT "# $cmd\n";
    }
  }

# prepare off insert vec fasta
  {
    my $in = $design->get_sample_off_insert_fasta($sample);
    my $list = $design->get_sample_off_insert_vec_inames($sample);
    my $out = $design->get_sample_off_insert_vec_fasta($sample);
    my $cmd = "$filter_fasta_script -i $BASE_DIR/$in -l $WORK_DIR/$list -o $BASE_DIR/$out";
		$cmd .= "\n$samtools faidx $BASE_DIR/$out";

    if(!-e "$BASE_DIR/$out") {
      print OUT "$cmd\n";
    }
    else {
			print STDERR "Warning: $BASE_DIR/$out already exists, won't override\n";
			$cmd =~ s/\n/\n# /sg;
      print OUT "# $cmd\n";
    }
  }

# prepare target insert ref map sorted and index
  {
    my $in = $design->get_sample_ref_map_filtered_sorted_file($sample);
    my $list = $design->get_sample_target_insert_vec_rnames($sample);
    my $out = $design->get_sample_target_insert_ref_sorted_file($sample);
		
    my $cmd = "if [ -s $WORK_DIR/$list ]; then java -jar $SCRIPT_DIR/$picard FilterSamReads --FILTER includeReadList -I $BASE_DIR/$in -O $BASE_DIR/$out --READ_LIST_FILE $WORK_DIR/$list;";
		$cmd .= "\nelse $samtools view -H $BASE_DIR/$in -b -o $BASE_DIR/$out; fi;";
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

# prepare off insert ref map sorted and index
  {
    my $in = $design->get_sample_ref_map_filtered_sorted_file($sample);
    my $list = $design->get_sample_off_insert_vec_rnames($sample);
    my $out = $design->get_sample_off_insert_ref_sorted_file($sample);

    my $cmd = "if [ -s $WORK_DIR/$list ]; then java -jar $SCRIPT_DIR/$picard FilterSamReads --FILTER includeReadList -I $BASE_DIR/$in -O $BASE_DIR/$out --READ_LIST_FILE $WORK_DIR/$list;";
		$cmd .= "\nelse $samtools view -H $BASE_DIR/$in -b -o $BASE_DIR/$out; fi;";
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

# prepare target genomic ref map sorted and index
  {
    my $in = $design->get_sample_ref_map_target_sorted_file($sample);
    my $list = $design->get_sample_target_insert_vec_rnames($sample);
    my $out = $design->get_sample_target_genomic_ref_sorted_file($sample);
		
    my $cmd = "if [ -s $WORK_DIR/$list ]; then java -jar $SCRIPT_DIR/$picard FilterSamReads --FILTER excludeReadList -I $BASE_DIR/$in -O $BASE_DIR/$out --READ_LIST_FILE $WORK_DIR/$list;";
		$cmd .= "\nelse $samtools view -h $BASE_DIR/$in -b -o $BASE_DIR/$out; fi;";
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

# prepare extract target genomic seq and index cmd
  {
    my $in = $design->get_sample_target_genomic_ref_sorted_file($sample);
    my $out = $design->get_sample_target_genomic_ref_map_seq($sample);
		
    my $cmd = "$samtools view -F 0x900 $BASE_DIR/$in -b | $samtools fasta -0 $BASE_DIR/$out $BASE_DIR/$in";
		$cmd .= "\n$samtools faidx $BASE_DIR/$out";

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
