#!/usr/bin/env perl
#########################################################################
# File Name: RNAfold.pl
# Author: Zhanghena
# mail: zhanghena1206@gmail.com
# Created Time: Wed 14 Mar 2018 09:41:22 PM CST
#########################################################################

use strict;
use warnings;
use Getopt::Long;

system "echo `date` dir: `pwd` command:'$0 @ARGV' >>$0.arg";
#args: perl RNAfold.pl *.bed base-pair probability

my $cutoff = shift;
my $BEDPath = shift;
my %pairs;
my @pair;


#ln -s /picb/rnomics3/Zhanghena/project/m5C/03_BSseq/03_merge_12_samples/03_noERCC_notRNA_BCE.bed
#perl -alne 'print join("\t",$F[0],$F[2],@F[2..5])' 03_noERCC_notRNA_BCE.bed >03_noERCC_notRNA_BCE.avinput
#annotate_variation.pl --geneanno --dbtype refGene --buildver hg38 --outfile 03_noERCC_notRNA_BCE.annoted 03_noERCC_notRNA_BCE.avinput humandb/
#cut -f 1,2 03_noERCC_notRNA_BCE.annoted.variant_function |grep -v intergenic|grep -v upstream|grep -v downstream|cut -f 2|perl -alne '$F[0] =~ s/\(.*//; @a=split(/,/,$F[0]); print join("\n",@a) ' |sort|uniq > 03_noERCC_notRNA_BCE.genes
#grep -f 03_noERCC_notRNA_BCE.genes humandb/hg38_refGene.txt > 03_noERCC_notRNA_BCE.transcripts
#cut -f 2- 03_noERCC_notRNA_BCE.transcripts | genePredToBed stdin 03_noERCC_notRNA_BCE.transcripts.bed
#bedtools getfasta -fi /picb/rnomics1/database/Human/hg38/genome/hg38_all.fa -bed 03_noERCC_notRNA_BCE.transcripts.bed -s -split -name >03_noERCC_notRNA_BCE.transcripts.fa

#split -l 26 ../03_noERCC_notRNA_BCE.transcripts.fa -d -a 3 list
#RNAfold -p --MEA=0.1 -T 75 --maxBPspan= -i 03_noERCC_notRNA_BCE.transcripts.fa > 03_noERCC_notRNA_BCE.RNAfold.out



#system('bedtools getfasta -fi /picb/rnomics1/database/Human/hg38/genome/hg38_all.fa -bed humandb/hg38_refGene.bed -s -split -name |head >test.fa');
#system('RNAfold -p --MEA=0.1 -T 75 --maxBPspan=150 -i test.fa > test.out');
#system('head -n 5 humandb/hg38_refGene.bed > test.bed');
my %BED = readBEDFile($BEDPath);
foreach my $KEY (keys %BED){
    my ($chr,$chromStart,$chromEnd,$name,$score,$strand,$thickStart,
        $thickEnd,$itemRgb,$blockCount,$blockSizes,$blockStarts) = split("\t",$BED{$KEY});
    my $dp = $name."__".$chr."_".$chromStart."-".$chromEnd."(".$strand.")"."_dp.ps";
    open DP, "< $dp" or die "can t open $dp\n";
    while(<DP>)
    {
        if (m/lbox/){
        @pair = split(/ /,$_);
        $pairs{$pair[0]} = $pair[1];
        }
    }
    close (DP);
    open DP, "< $dp" or die "can't open $dp\n";
    my @sites = ();
    while(<DP>)
    {
        if (m/ubox/){
        @pair = split(/ /,$_);
            if (exists $pairs{$pair[0]}){
                if ($pair[1] == $pairs{$pair[0]} && $pair[2] >= $cutoff){
                    push(@sites, $pair[0]);
                    push(@sites, $pair[1]);
                }
            }
        }
    }
    close (DP);
    my @genomic = ();
    foreach my $i (@sites){
        my $site = transcript_to_genomic($i,$BED{$KEY});
        push(@genomic, $site);
}
    
}

########################################################################## transform the transcript position to genomic position
sub transcript_to_genomic{
    #0-based
    #if strand is "-", return position of strand.
    #---my $genomeic = transcript_to_genomic()
    my $transcript_pos = shift;
    my $BED = shift;
    print $BED;
    my ($chr,$chromStart,$chromEnd,$name,$score,$strand,$thickStart,
        $thickEnd,$itemRgb,$blockCount,$blockSizes,$blockStarts) = split("\t",$BED);
    my $KEY = ">".$name."::".$chr.":".$chromStart."-".$chromEnd."(".$strand.")" ;
    my @blockSize   = split(/,/,$blockSizes);
    my @blockStart  = split(/,/,$blockStarts);
    my $length      = 0;
    my $genomic_pos = 0;
    if ($strand eq "+"){
        foreach my $i (0..$blockCount-1){
            if ($transcript_pos < $length + $blockSize[$i]){
                $genomic_pos = $chromStart + $blockStart[$i] + $transcript_pos - $length ; 
                print $chromStart."\t".$blockStart[$i]."\t".$transcript_pos."\t".$KEY;
                last;
            }
            $length += $blockSize[$i]; 
        }
    }
    else {
        for(my $i = $blockCount-1;$i>=0;$i--){
            if ($transcript_pos < $length + $blockSize[$i]){
                $genomic_pos = $chromStart + $blockStart[$i] + $blockSize[$i] - ($transcript_pos - $length);
                last;
            }
            $length += $blockSize[$i]; 
        }
    }
    print $chr."\t".$genomic_pos."\t".$strand."\t".$KEY."\t".$transcript_pos."\n";
}


########################################################################## readBEDFile
sub readBEDFile {
    #---my $BEDInfoHsh_ref = readBEDFile($BEDPath);
    my $BEDPath = $_[0];
    my %BEDHash;
    open BED, "< $BEDPath" or die "can t open $BEDPath\n";
    while (my $theLine = <BED>) {
        chomp $theLine;
        #chr19    47045908    47048624    NM_017854    0    -    47045986    47048614    0    3    344,93,218,    0,684,2498,
        my ($chr,$chromStart,$chromEnd,$name,$score,$strand,$thickStart,
        $thickEnd,$itemRgb,$blockCount,$blockSizes,$blockStarts) = split("\t",$theLine);
        my $KEY = ">".$name."::".$chr.":".$chromStart."-".$chromEnd."(".$strand.")" ;
        $BEDHash{$KEY} = $theLine;
        
    }
    close BED;
    return %BEDHash;
}
