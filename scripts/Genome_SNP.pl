#!/usr/bin/perl
use strict;
use warnings;

my $ref=$ARGV[0];#
my $f=$ARGV[1];#Genome.list

open(FILE,$f);my %snp;my %site;my %ref;
while(1){
	my $l=<FILE>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	unless(-e $a[1]){next;}
	my $type=0;my %high;
	if(-e $a[2]){
		$type=1;
		open(S,$a[2]);
		
		while(1){
			my $line=<S>;
			unless($line){last;}
			chomp $line;
			my @b=split"\t",$line;
		#	unless($b[0] eq "NC_018042.1"){next;}
			my @c=split";",$b[7];
			my $dp=substr($c[0],3);
			
			if(length($b[3])==1 && length($b[4])==1 && $b[5]>=50 && $dp>=10){
				my $site=$b[0].",".$b[1];
				$high{$site}=$b[4];$snp{$a[0]}{$site}=$b[2];$site{$site}++;$ref{$site}=$b[1];
				
			}
		}
		close S;
	}
	open(M,$a[1]);
		while(1){
			my $line=<M>;
			unless($line){last;}
			chomp $line;
			my @d=split"\t",$line;
		#	unless($d[10]=~/NC_018042.1/){next;}
			unless($d[1]=~/A|T|G|C/ && $d[2]=~/A|T|G|C/){next;}
			my $site=$d[10].",".$d[0];
			if($type==1){
				if(exists $high{$site}){$snp{$a[0]}{$site}=$d[2];$site{$site}++;$ref{$site}=$d[1];}
			}else{
				$snp{$a[0]}{$site}=$d[2];$site{$site}++;$ref{$site}=$d[1];
			}
		}
		close M;
}
close FILE;

open(OUT,">SNP_cat.list");my %seq;
my @sample=keys %snp;
my @site=sort keys %site;
print OUT "ID\tRef\t";
foreach my $sample (@sample) {
	print OUT "$sample\t";
}
print OUT "\n";
foreach my $site (@site) {
	print OUT "$site\t$ref{$site}\t";
	$seq{"Ref"}.=$ref{$site};
	foreach my $sample (@sample) {
		if(exists $snp{$sample}{$site}){print OUT "$snp{$sample}{$site}\t";$seq{$sample}.=$snp{$sample}{$site};}else{print OUT "$ref{$site}\t";$seq{$sample}.=$ref{$site};}
	}
	print OUT "\n";
}
close OUT;

open(OUT,">Genome_SNP.fasta");
my @seq=keys %seq;
foreach my $seq (@seq) {
	print OUT ">$seq\n$seq{$seq}\n";
}