#!/bin/env perl
# get per-insert annotation summary
use strict;
use warnings;

use Getopt::Long;

my $min_ratio = 0.9;
my $feat_tag = 'label';
my $usage = "Usage: $0 GFF-INFILE TSV-OUTFILE <--feat-basic STR> <--feat-full STR> [--min-ratio $min_ratio] [--feature-tag $feat_tag]";

# get opts
my $infile = shift or die $usage;
my $outfile = shift or die $usage;

my $feat_basic;
my $feat_full;

GetOptions(
"feat-basic=s" => \$feat_basic,
"feat-full=s" => \$feat_full,
"min-ratio=f" => \$min_ratio,
"feat-tag=s" => \$feat_tag)
or die "Error in command line arguments, usage: $usage";

defined($feat_basic) && defined($feat_full) && $min_ratio > 0 || die $usage;

# open files
open(IN, "<$infile") || die "Unable to open $infile: $!";
open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";

my @feat_basic = split(/\|/, $feat_basic);
my @feat_full = split(/\|/, $feat_full);

# read in GFF annotation
my %name2feat_cover;
while(my $line = <IN>) {
  chomp $line;
	next if($line =~ /^#/);
	my ($name, $src, $type, $start, $end, $score, $strand, $frame, $attrs) = split(/\t/, $line);
	my ($label, $cover_ratio);
	foreach my $attr (split(/;/, $attrs)) {
		my ($tag, $val) = split(/=/, $attr);
		if($tag eq $feat_tag) {
			$label = $val;
		}
		elsif($tag eq 'CoverRatio') {
			$cover_ratio = $val;
		}
	}
	push(@{$name2feat_cover{$name}}, [$label, $cover_ratio]);
}

# scan each insert and output
# print OUT "#options invoked: --feat-basic=$feat_basic;--feat-full=$feat_full;--min-ratio=$min_ratio\n";
print OUT "insert_id\tfeature_cover_summ\tdonor_vec_basic_clone\tdonor_vec_full_clone\n";

foreach my $insert_id (sort keys %name2feat_cover) {
  my $cover_summ = join('|', map { join(':', @$_) } @{$name2feat_cover{$insert_id}});
  # filter feature with good cover-ratio
  my @feat_fwd;
  foreach my $feat (@{$name2feat_cover{$insert_id}}) {
    if($feat->[1] >= $min_ratio) {
      push(@feat_fwd, $feat->[0]);
    }
  }
	my @feat_rev = reverse(@feat_fwd);
# count feature fwd
  my $basic_clone = count_sub_array(\@feat_fwd, \@feat_basic);
  my $full_clone = count_sub_array(\@feat_fwd, \@feat_full);
  if(@feat_basic > 1) { # also count feature rev
    $basic_clone += count_sub_array(\@feat_rev, \@feat_basic);
  }
  if(@feat_full > 1) { # also count feature rev
    $full_clone += count_sub_array(\@feat_rev, \@feat_full);
  }
	print OUT "$insert_id\t$cover_summ\t$basic_clone\t$full_clone\n";
}

close(IN);
close(OUT);

sub count_sub_array {
  my ($A, $B) = @_;
  my $n = @$A;
  my $m = @$B;
  my $clone = 0;
  my ($i, $j) = (0, 0);
  while($i < $n && $j < $m) {
    if($A->[$i] eq $B->[$j]) {
      $i++;
      $j++;

      if($j == $m) {
        $clone++;
        $i++;
        $j = 0;
      }
    }
    else {
      $i = $i - $j + 1;
      $j = 0;
    }
  }

  return $clone;
}
