#!/usr/bin/env perl
# get various statistics for Cas9_LongRead runs
use strict;
use warnings;
# use File::Basename;

my $usage = "Usage: $0 MANIFEST-FILE OUTFILE";
my $manfile = shift or die $usage;
my $outfile = shift or die $usage;

my @header = qw(run_name study_ID total_read chr1_mapped vec_mapped chr1_and_vec_mapped total_ROI_read filtered_ROI_read donor_mapped_ROI_read nuclease_mapped_ROI_read donor_and_nuclease_mapped_ROI_read);

open(MAN, "<$manfile") || die "Unable to read $manfile: $!";
open(OUT, ">$outfile") || die "Unable to write to $outfile: $!";
print OUT join("\t", @header), "\n";

# process each run
<MAN>;
while(my $line = <MAN>) {
	chomp $line;
	my ($runname, $studyID, $nuclease_vec, $donor_vec) = split(/\t/, $line);
	print STDERR "Getting stats for $runname ($studyID)\n";
# get total read
	my $total_read;
	{
		my $in = "$studyID\_all_pass.fastq";
		$total_read = `cat $run/$in | wc -l`;
		chomp $total_read;
		$total_read /= 4;
	}

# get chr1_mapped, nuclease_mapped and donor mapped
	my $chr1_mapped;
	{

	}

  print OUT "$runname\t$studyID\t$total_read\n";
}

close(OUT);
