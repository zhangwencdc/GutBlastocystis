#!/usr/bin/perl
use strict;
#use warnings;

my $ref=$ARGV[0];#
my $f=$ARGV[1];#CHMP1.Mito.list
###与virus的比对参数 可调整
my $cutoff;my $match_cutoff;
if ( !defined $match_cutoff ){ $match_cutoff=70;}###最小比对长度70
if ( !defined $cutoff ){ $cutoff=80; }##最小比对百分百80%
open(FILE,$f);my %snp;my %site;my %ref;my $num;my %name;$name{"Ref"}="Ref";
open(OH,">Candidate_Hetero_SNP.list");
while(1){
	my $l=<FILE>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	#unless(-e $a[1]){next;}
	my $type=0;
	$name{$a[0]}=$a[1];
		$num++;
		$type=1;
		my $fq1=$a[2];my $fq2=$a[3];
	#	print "/home/zhangwen/bin/bowtie2-2.4.4-linux-x86_64/bowtie2 -1 $fq1 -2 $fq2 -x $ref -S $a[0].target.sam --no-unal\n";

		###vcf
		open(S,"$a[0].target.vcf")||print "Error:$l\n";
		my $he;
		while(1){
			my $line=<S>;
			unless($line){last;}
			chomp $line;
			my @b=split"\t",$line;
			#unless($b[0]=~/JZRK/){next;}
			my @c=split";",$b[7];
			my $dp=substr($c[0],3);
			unless ($b[4]=~/A|T|G|C/) {
				$b[4]=$b[3]; #将. 替换成ref site
			}
			if(length($b[3])==1 && length($b[4])==1 && $b[5]>=50 && $dp>=5){
				my $site=$b[0].",".$b[1];
				$snp{$a[0]}{$site}=$b[4];$site{$site}++;$ref{$site}=$b[3];
				my $t=substr($b[9],0,3);
				if($t eq "0/1" || $t eq "1/0"){
				if($b[3] ne $b[4]){print OH "$a[0] $line\n";$he++;}
				}
			}
		}
		close S;print "Canddiate Heter SNP,$a[0],$he\n";
	
	
}
close FILE;
print"Sample Num: $num\n";
open(OUT,">SNP_Mito_cat.list");my %seq;
my @sample=keys %snp;
my @site=sort keys %site;
print OUT "ID\tRef\t";
foreach my $sample (@sample) {
	print OUT "$name{$sample}\t";
}
print OUT "\n";
foreach my $site (@site) {
	print OUT "$site\t$ref{$site}\t";
	$seq{"Ref"}.=$ref{$site};my $mut;my %type;
	foreach my $sample (@sample) {
		if(exists $snp{$sample}{$site}){print OUT "$snp{$sample}{$site}\t";$seq{$sample}.=$snp{$sample}{$site};
		if($snp{$sample}{$site} ne $ref{$site}){$mut++;}$type{$snp{$sample}{$site}}++;
		}else{print OUT "-\t";$seq{$sample}.="-";}
	}
	if($mut>0.2*$num && $mut<($num-0.2*$num)){print OUT "Candidate_for_evolution,";}
	my @type=keys %type;my $typen=@type;if($typen>=2){print OUT "Within_mut,";}
	print OUT "\n";
}
close OUT;

open(OUT,">Mito_SNP.fasta");
my @seq=keys %seq;
foreach my $seq (@seq) {
	print OUT ">$name{$seq}\n$seq{$seq}\n";
}