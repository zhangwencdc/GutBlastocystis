#!/usr/bin/perl
use strict;
use warnings;
#Usage: 基于Meta Assemble结果，得到Blastocystis的线粒体序列
use FindBin qw($Bin $Script);
use Getopt::Long;

my ($file,$outdir,$mark,$Ref,$Ref_s,$identity,$coverage,$type,$HELP,$NUMREPS,$BURNIN,$sum_n,$STEP,$align_cut,$key);
GetOptions(
                "I:s"=>\$file,  # Meta Assembled Genome
                "O:s"=>\$outdir,     ###outdir
                "K|Key|key|tag:s"=>\$key,
                "R|Ref:s"=>\$Ref,  #Blastocystis_sub3_Mitochondria.fasta   
                "A|Align:i"=>\$align_cut, 
                "C|Coverage:i"=>\$coverage, 
            
                "help"=>\$HELP
);
die `pod2text $0` if ($HELP || !defined $file );

if(!defined $Ref){$Ref=$Bin."/Blastocystis-sub3_Genome.fasta";}
#if(!defined $coverage){$coverage=0.5;}  #默认Coverage cutoff 50%
#if(!defined $align_cut){$align_cut=20000;} #默认 最少algin 20kb
if(!defined $key){$key="input";}
if(!defined $outdir){$outdir="./";}else{system "mkdir $outdir\n";}
###
my $dnadiff="dnadiff";
my $seqkit="seqkit";
print "$dnadiff $Ref $file -p tmp\n";

system "$dnadiff $Ref $file -p tmp\n";

my $cf="tmp.1coords";
	my $align;
	open(C,$cf);
	while(1){
		my $l=<C>;
		unless($l){last;}
		chomp $l;
		my @b=split"\t",$l;
		unless($b[4]>=1000){next;}  #仅考虑align>1000bp的比对结果
		system "$seqkit grep -p $b[12] $file >>tmp.fas\n";
		$align+=$b[4];
	
	}
	close C;
	print "$key,$align bp\n";
	system "$seqkit rmdup tmp.fas >$outdir/$key.target.fas\n";
	system "rm -rf tmp*\n";
	system "$dnadiff $Ref $outdir/$key.target.fas -p $outdir/$key\n";