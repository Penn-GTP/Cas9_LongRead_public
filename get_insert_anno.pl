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

my $options = qq([OPTIONS]
OPTIONS:
  --include-type [TYPE1[,TYPE2]]: include given types in the GFF file as the annotation source, multiple values allowed separated by comma or by given multiple times
  --exclude-type [TYPE3[,TYPE4]]: exclude given types in the GFF file as the annotation source, multiple values allowed separated by comma or by given multiple times; [default: ) . join(",", @exc_types) . qq(]
  --display [TAG]: use TAG as the displayName tag for better 3rd party tool (i.e. IGV) visulization; [default: $display]
);

my $usage = "Usage: $0 FEATURE-GFFFILE BED-INFILE GFF-OUTFILE $options";

my $gff_file = shift or die $usage;
my $infile = shift or die $usage;
my $outfile = shift or die $usage;

GetOptions(
"include-type=s" => \@inc_types,
"exclude-type=s" => \@exc_types,
"display=s" => \$display)
or die "Error in command line arguments, usage: $usage";

# expand potential comma-separated opts
@inc_types = split(/,/, join(',', @inc_types));
@exc_types = split(/,/, join(',', @exc_types));

open(IN, "bedtools intersect -a $infile -b $gff_file -wo -bed |") || die "Unable to open $infile: $!";
open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";

# read in insert anno
my %insert2overlap;

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
		my ($q_len) = $name =~ /:(\d+)[IS]:/;

# build target2query index, all coordinates are 0-based
		my @r2q_idx = 0 x $t_len;
		my $i = 0; # index on target
		my $j = 0; # index on query
		while($cigar =~ /(\d+)([MIDNSHPX=])/g) {
			my ($len, $op) = ($1, $2);
			for(my $k = 0; $k < $len; $k++) {
				$r2q_idx[$i] = $j;
				if($op eq 'M' || $op eq '=' || $op eq 'X' || $op eq 'D' || $op eq 'N') {
					$i++;
				}
				if($op eq 'M' || $op eq '=' || $op eq 'X' || $op eq 'I' || $op eq 'S') {
					$j++;
				}
			}
		}

		my $over_start = $start > $start2 ? $start : $start2;
		my $over_end = $end < $end2 ? $end : $end2;

		my $over_from = $r2q_idx[$over_start - $start] + 1; # 1-based for GFF output
		my $over_to = $r2q_idx[$over_end - 1 - $start] + 1;

		my $over_strand = $strand eq $strand2 ? '+' : '-';

		if($strand eq '-') { # query is reverse complemented
			($over_from, $over_to) = ($q_len - $over_to + 1, $q_len - $over_from + 1);
		}

		my $over_attr = $attr;

# add overlap info
		push(@{$insert2overlap{$name}}, [$name, $src, $type, $over_from, $over_to, $score, $over_strand, $frame, $over_attr]);
	}
}

# output
print OUT "##gff-version 3\n";
print OUT "##displayName=$display\n";

foreach my $name (sort keys %insert2overlap) {
	my %id2ver;
  foreach my $overlap (sort { $a->[3] <=> $b->[3] } @{$insert2overlap{$name}}) {
    # only show overlap features with IDs available and not a whole chromosome
		if($overlap->[8] =~ /ID=([^;]+)/) {
			my $id = !exists $id2ver{$1} ? $1 : $1 . '.' . $id2ver{$1};
			$id2ver{$1}++;
			$overlap->[8] =~ s/ID=([^;]+)/ID=$id/;
			print OUT join("\t", @$overlap), "\n";
		}
	}
}

close(IN);
close(OUT);
