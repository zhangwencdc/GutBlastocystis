#!/usr/bin/perl
use strict;
use warnings;

my $ref=$ARGV[0];#
my $f=$ARGV[1];#Mito.list

open(FILE,$f);my %snp;my %site;my %ref;
open(OH,">Candidate_Hetero_SNP.list");
while(1){
	my $l=<FILE>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	print "$a[2]\n";
	unless(substr($a[2],length($a[2])-1,1)=~/[0-9a-zA-Z]/){$a[2]=substr($a[2],0,length($a[2])-1);}
	#unless(-e $a[1]){next;}
	my $type=0;my %high;
	if(-e $a[2]){
		$type=1;
		open(S,$a[2]);
		open (SL,">$a[0].SNP.list");my $he=0;
		while(1){
			my $line=<S>;
			unless($line){last;}
			chomp $line;
			my @b=split"\t",$line;
			unless($b[0] eq "NC_018042.1"){next;}
			my @c=split";",$b[7];
			my $dp=substr($c[0],3);
			
			if(length($b[3])==1 && length($b[4])==1 && $b[5]>=50 && $dp>=10){
				$high{$b[1]}=$b[4];my $tmp=$b[1]+1;
				#print "Mitochondrion\t$b[1]\t$tmp\tfill_color=vlpblack\n";
				print SL "Mitochondrion\t$b[1]\t$tmp\tfill_color=vlpblack\n";
				my $t=substr($b[9],0,3);
					if($t eq "0/1" || $t eq "1/0"){  ###潜在的杂合位点
					print "$line\n";
				print OH "$a[0] $line\n";$he++;
					}
			}
		}
		close S;close SL;print "Canddiate Heter SNP,$a[0],$he\n";
	}
	open(M,$a[1]);
		while(1){
			my $line=<M>;
			unless($line){last;}
			chomp $line;
			my @d=split"\t",$line;
			unless($d[10]=~/NC_018042.1/){next;}
			unless($d[1]=~/A|T|G|C/ && $d[2]=~/A|T|G|C/){next;}
			if($type==1){
				if(exists $high{$d[0]}){$snp{$a[0]}{$d[0]}=$d[2];$site{$d[0]}++;$ref{$d[0]}=$d[1];}
			}else{
				$snp{$a[0]}{$d[0]}=$d[2];$site{$d[0]}++;$ref{$d[0]}=$d[1];
			}
		}
		close M;
}
close FILE;

open(OUT,">SNP_cat.list");my %seq;
my @sample=keys %snp;
my @site=sort{$a<=>$b} keys %site;
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

open(OUT,">Mito_SNP.fasta");
my @seq=keys %seq;
foreach my $seq (@seq) {
	print OUT ">$seq\n$seq{$seq}\n";
}