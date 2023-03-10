#!/usr/bin/env perl
# This script is ued to get insert annotation details from anno BED files
use strict;
use warnings;
use sort 'stable';
use Getopt::Long;
use List::MoreUtils qw(any none);

my @inc_types;
my @exc_types = qw(source chromosome region);
my $display = 'Name';

my $exc_types = join(",", @exc_types);
my $options = qq([OPTIONS]
OPTIONS:
  --include-type [TYPE1[,TYPE2]]: include given types in the GFF file as the annotation source, multiple values allowed separated by comma or by given multiple times
  --exclude-type [TYPE3[,TYPE4]]: exclude given types in the GFF file as the annotation source, multiple values allowed separated by comma or by given multiple times; [default: $exc_types]
  --display [TAG]: use TAG as the displayName tag for better 3rd party tool (i.e. IGV) visulization; [default: $display]
  --genome-size FILE: genome size file of 2 columns, with 1st field chromosome names and 2nd field their sizes [default: null]);

my $usage = "Usage: $0 GFF-INFILE BED-INFILE OUTFILE $options";

my $gfffile = shift or die $usage;
my $bedfile = shift or die $usage;
my $outfile = shift or die $usage;
my $genome_file;

GetOptions(
"include-type=s" => \@inc_types,
"exclude-type=s" => \@exc_types,
"display=s" => \$display,
"genome-size=s" => \$genome_file)
or die "Error in command line arguments, usage: $usage";

# expand potential comma-separated opts
@inc_types = split(/,/, join(',', @inc_types));
@exc_types = split(/,/, join(',', @exc_types));
my %name2overlap;
my $overlap_id = '';
my @r2q_idx;
my @r2op_idx;

my $cmd = "bedtools intersect -a $bedfile -b $gfffile -wo -bed";
if($genome_file) {
	$cmd .= "-g $genome_file";
}

open(IN, "$cmd |") || die "Unable to open $bedfile: $!";
open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";

# read in insert anno
while(my $line = <IN>) {
	chomp $line;
	my @fields = split(/\t/, $line);
	if(@fields < 17) {
		print STDERR "input BED file must contains the 7th column of the cigar string of the alignment\n";
		exit;
	}
	my ($chr, $start, $end, $name, $score, $strand, $cigar, $chr2, $src, $type, $start2, $end2, $score2, $strand2, $frame, $attr, $over_len) = @fields;
	$start2--;
	if((any { $type eq $_ } @inc_types) || (none { $type eq $_ } @exc_types) ) {
		my $t_len = $end - $start;
		my ($q_len) = get_qlen_by_cigar($cigar);

# build target2query relative-pos index, all coordinates are 0-based
    if($overlap_id ne "$chr:$start:$end:$name:$cigar") {
			@r2q_idx = (0) x $t_len;
			@r2op_idx = (0) x $t_len;
			my $i = 0;
			my $j = 0;
			while($cigar =~ /(\d+)([MIDNSHPX=])/g) {
				my ($len, $op) = ($1, $2);
				for(my $k = 0; $k < $len; $k++) {
					$r2q_idx[$i] = $j;
					$r2op_idx[$i] = $op;
					if($op eq 'M' || $op eq '=' || $op eq 'X' || $op eq 'D' || $op eq 'N') {
						$i++;
					}
					if($op eq 'M' || $op eq '=' || $op eq 'X' || $op eq 'I' || $op eq 'S') {
						$j++;
					}
				}
			}
		}
		$overlap_id = "$chr:$start:$end:$name:$cigar";

		my $over_start = $start > $start2 ? $start : $start2;
		my $over_end = $end < $end2 ? $end : $end2;

		my $over_from = $r2q_idx[$over_start - $start] + 1;
		my $over_to = $r2q_idx[$over_end - 1 - $start] + 1;

		my $over_strand = $strand eq $strand2 ? '+' : '-';

		if($strand eq '-') {
			($over_from, $over_to) = ($q_len - $over_to + 1, $q_len - $over_from + 1);
		}

		# calculate overlap bases
		my $over_base = 0;
		for(my $i = $over_start - $start; $i < $over_end - $start; $i++) {
			if($r2op_idx[$i] eq 'M' || $r2op_idx[$i] eq '=' || $r2op_idx[$i] eq 'X') {
				$over_base++;
			}
		}
    my $over_ratio = $over_base / ($end2 - $start2 + 1);

# build overlap attr
    my $over_attr = "$attr;CoverLen=$over_len;CoverBase=$over_base;CoverRatio=$over_ratio";

# add overlap info
		push(@{$name2overlap{$name}}, [$name, $src, $type, $over_from, $over_to, $score, $over_strand, $frame, $over_attr]);
	}
}

# output
print OUT "##gff-version 3\n";
print OUT "##displayName=$display\n";

foreach my $name (sort keys %name2overlap) {
	my %id2ver;
	foreach my $overlap (sort { $a->[3] <=> $b->[3] } @{$name2overlap{$name}}) {
# make ID unique in features
		if($overlap->[8] =~ /ID=([^;]+)/) {
			my $id = !exists $id2ver{$1} ? $1 : $1 . '.' . $id2ver{$1};
			$id2ver{$1}++;
			$overlap->[8] =~ s/ID=([^;]+)/ID=$id/;
		}
		print OUT join("\t", @$overlap), "\n";
	}
}

close(IN);
close(OUT);

sub get_qlen_by_cigar {
	my $cigar = shift;
	my $qlen = 0;
	while($cigar =~ /(\d+)([MIDNSHPX=])/g) {
		my ($len, $op) = ($1, $2);
		if($op eq 'M' || $op eq 'I' || $op eq 'S' || $op eq '=' || $op eq 'X') {
			$qlen += $len;
		}
	}
	return $qlen;
}
