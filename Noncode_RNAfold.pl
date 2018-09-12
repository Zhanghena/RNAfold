#!/usr/bin/env perl
#########################################################################
# File Name: Noncode_RNAfold.pl
# Author: Zhenghena
# mail: zhanghena@gmail.com
# Created Time: Thu 09 Aug 2018 10:34:44 AM CST
#########################################################################

use strict;
use warnings;
use Getopt::Long;

sub usage {
	print <<"END_USAGE";
Usage: perl $0
	--input tab fasta file of 300 nt around m5C sites
	--output pairs sites
	--p probability(default >=0.8)
END_USAGE
	exit;
}

my ($input,$output,$p);
$p=0.8;

GetOptions (
	'input=s' =>\$input,
	'output=s'=>\$output,
	'p=f'     =>\$p
) or usage();
usage() if (!$input or !$output);

open my $in,   "< $input" or die;
##chr19   47048329        47048630        0       0       -		TGGGGCGGTACGGCGCGCGGCGGCCGTAGTGCGAATCATGGGAGGCGGCTGGTGGTGGGCTCGGGCCGCTCGCCTTGCCCGTCTTCGCTTCCGGAGGTCGCTACTGCCGCCTCAGCGGCCCCGGAGCGGGGGCGCCCGGGGGTCCTTCGCCCCCGGCCACGGTCCCCGCGCCGGGGCTTCGCCGCCCCCAGTGTCCGAGCTGGATCGTGCGGACGCCTGGCTCCTCCGAAAAGCGCACGAGACAGGTCGGTGGCTAGGGGTGGGGGCGGTGCGGATGGGACGGGGGTCCCGGGCTCGCGTG

open my $out,  "> $output" or die;

while(<$in>){
	my @F=split /\t/;
	my $loci=$F[0].":".($F[2]-150).":".$F[5];
	my $seq=$F[6];
	
	open my $temp, "> temp.fa" or die;
	print $temp ">".$loci."\n".$seq;
	close $temp;
	
	system("RNAfold -p --MEA=0.1 -T 75 --maxBPspan=150 -i temp.fa --id-prefix temp >temp.log");

	open (my $DP, "< temp_0001_dp.ps") or die "can not open file\n";
	my %pairs;
	while (<$DP>){
		if (m/lbox/ && m/151/ ) {
			my @pair = split(/ /,$_);
			$pairs{$pair[0]} = $pair[1];
		}
	}
	close $DP;

	open ($DP, "< temp_0001_dp.ps") or die "can not open file\n";
	while(<$DP>){
		if (m/ubox/){
			my @pair = split(/ /,$_);
			if (exists $pairs{$pair[0]}){
				if ($pair[1] == $pairs{$pair[0]} && $pair[2] >= $p){
					print $out $loci."\t".$pair[2]."\n";
				}	
			}
		}
	}
	close $DP;

}
close $in;
close $out;







