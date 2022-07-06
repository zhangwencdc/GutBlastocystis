#!/usr/bin/perl
use strict;
#use warnings;

my $ref=$ARGV[0];#
my $f=$ARGV[1];#CHMP1.list

open(FILE,$f);my %snp;my %site;my %ref;my $num;my %name;$name{"Ref"}="Ref";
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
		open(S,$a[2])||print "Error:$l\n";
		
		while(1){
			my $line=<S>;
			unless($line){last;}
			chomp $line;
			my @b=split"\t",$line;
			unless($b[0]=~/JZRK/){next;}
			my @c=split";",$b[7];
			my $dp=substr($c[0],3);
			
			if(length($b[3])==1 && length($b[4])==1 && $b[5]>=50 && $dp>=10){
				my $site=$b[0].",".$b[1];
				$snp{$a[0]}{$site}=$b[4];$site{$site}++;$ref{$site}=$b[3];
				
			}
		}
		close S;
	
	
}
close FILE;
print"Sample Num: $num\n";
open(OUT,">SNP_cat.list");my %seq;
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
		if(exists $snp{$sample}{$site}){print OUT "$snp{$sample}{$site}\t";$seq{$sample}.=$snp{$sample}{$site};$mut++;$type{$snp{$sample}{$site}}++;}else{print OUT "$ref{$site}\t";$seq{$sample}.=$ref{$site};$type{$ref{$site}}++;}
	}
	if($mut>0.2*$num && $mut<($num-0.2*$num)){print OUT "Candidate_for_evolution,";}
	my @type=keys %type;my $typen=@type;if($typen>2){print OUT "Tri_mut,";}
	print OUT "\n";
}
close OUT;

open(OUT,">Genome_SNP.fasta");
my @seq=keys %seq;
foreach my $seq (@seq) {
	print OUT ">$name{$seq}\n$seq{$seq}\n";
}