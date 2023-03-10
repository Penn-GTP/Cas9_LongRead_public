#!/bin/env perl
# get per-insert annotation summary
use strict;
use warnings;

use Getopt::Long;

my $min_ratio = 0.9;
my $feat_tag = 'label';
my $ARM_key = "(?i:ARM)|shHDR";
my $ARM_min_ratio = 0.1;

my $usage = "Usage: $0 GFF-INFILE TSV-OUTFILE <--func-feat STR> [--min-ratio $min_ratio] [--feature-tag $feat_tag] [--ARM-key $ARM_key] [--ARM-min-ratio $ARM_min_ratio]";

# get opts
my $infile = shift or die $usage;
my $outfile = shift or die $usage;

my $func_feat;

GetOptions(
"func-feat=s" => \$func_feat,
"min-ratio=f" => \$min_ratio,
"feat-tag=s" => \$feat_tag,
"ARM-key=s" => \$ARM_key,
"ARM-min-ratio=f" => \$ARM_min_ratio)
or die "Error in command line arguments, usage: $usage";

defined($func_feat && $min_ratio > 0 && $ARM_min_ratio > 0) || die $usage;

# open files
open(IN, "<$infile") || die "Unable to open $infile: $!";
open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";

my @func_feat = split(/\|/, $func_feat);
$ARM_key = qr/$ARM_key/;

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
	if(defined $label) {
		push(@{$name2feat_cover{$name}}, [$label, $cover_ratio]);
	}
}

# scan each insert and output
# print OUT "#options invoked: --feat-basic=$feat_basic;--feat-full=$feat_full;--min-ratio=$min_ratio\n";
print OUT "rname\tfeature_cover_summ\tfunctional_clone\tinsert_type\n";

foreach my $rname (sort keys %name2feat_cover) {
	my @all_feats = @{$name2feat_cover{$rname}};
  my $cover_summ = join('|', map { join(':', @$_) } @all_feats);
  # filter feature with good cover-ratio
  my @feat_fwd;
  foreach my $feat (@all_feats) {
    if($feat->[1] >= $min_ratio) {
      push(@feat_fwd, $feat->[0]);
    }
  }
	my @feat_rev = reverse(@feat_fwd);
# count feature fwd
  my $func_clone = count_sub_array(\@feat_fwd, \@func_feat);
  if(@func_feat > 1) { # also count feature rev
    $func_clone += count_sub_array(\@feat_rev, \@func_feat);
  }
  my $insert_type = ends_with_ARM($ARM_key, $ARM_min_ratio, @all_feats) ? 'HDR' : 'NHEJ';
	print OUT "$rname\t$cover_summ\t$func_clone\t$insert_type\n";
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

sub ends_with_ARM {
  my ($ARM_key, $ARM_min_ratio, @all_feats) = @_;
  return @all_feats > 0 && ($all_feats[0][0] =~ /$ARM_key/ && $all_feats[0][1] >= $ARM_min_ratio || $all_feats[$#all_feats][0] =~ /$ARM_key/ && $all_feats[$#all_feats][1] >= $ARM_min_ratio);
}

sub contains_ARM {
  my ($ARM_key, @all_feats) = @_;
  return scalar grep { $_->[0] =~ /$ARM_key/ } @all_feats;
}
