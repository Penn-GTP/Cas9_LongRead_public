#!/bin/env perl
use strict;
use warnings;
# This script is used to filter fasta file given a filter geneID list

my $usage = "Usage: $0 -i FASTA-INFILE -l FILTER-ID-LIST -o OUTFILE [-v(inverse)] [--id-regex (\\S+)]";
my $infile;
my $filter_file;
my $outfile;
my $inverse = 0;
my $id_regex = qr/(\S+)/;

for(my $i = 0; $i < @ARGV; $i++) {
  if($ARGV[$i] eq '-i') {
    $infile = $ARGV[++$i];
  }
  elsif($ARGV[$i] eq '-l') {
    $filter_file = $ARGV[++$i];
  }
  elsif($ARGV[$i] eq '-o') {
    $outfile = $ARGV[++$i];
  }
  elsif($ARGV[$i] eq '-v') {
    $inverse = 1;
  }
  elsif($ARGV[$i] eq '--id-regex') {
    my $regex_str = $ARGV[++$i];
    $regex_str =~ s/^["']//; # remove leading quotes
    $regex_str =~ s/["']$//; # remove tailing quotes
    $id_regex = qr/$regex_str/;
    print STDERR "Using id regex: $id_regex\n";
  }
  else {
    print STDERR "Unknown option '$ARGV[$i]'.\n";
    exit;
  }
}

# Check options
unless(defined $infile && defined $filter_file && defined $outfile) {
  print STDERR "$usage\n";
  exit;
}
open(IN, "<$infile") || die "Unable to open $infile: $!";
open(FILTER, "<$filter_file") || die "Unable to open $filter_file: $!";
open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";

# Get filter name into a hash!
print STDERR "Read in filter list ...\n";
my %filter;
while(my $line = <FILTER>) {
  chomp $line;
  $filter{$line} = 1;
}

# Filter input file
print STDERR "Searching input file ...\n";
my $flag = 0;
while(my $line = <IN>) {
  if($line =~ /^>$id_regex/) {  # header line
    $flag = exists $filter{$1};
    if($inverse) {  # Filterout instead of filter
      $flag = !$flag;
    }
  }
  if($flag) {
    print OUT $line;
  }
}

print "Done!\n";
close(IN);
close(FILTER);
close(OUT);
