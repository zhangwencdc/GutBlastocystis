#!/usr/bin/perl
use strict;
#use warnings;
#Usage：SNP list for circos; SNP 同义突变、非同义突变数量；发生在gene上和gene间区的SNP数量统计
my $file=$ARGV[0];#/home/zhangwen/project/2022Time/WGS_Analysis/Eukaryote/New_Blastocystis/Blastocystis_Mitochondria 输入文件路径
my $cds=$ARGV[1];#Blastocystis_sub3_Mitochondria_cds-v2.fasta
my $out=$ARGV[2];#Mitochondria.SNP.stat
###读取cds信息
open(CDS,$cds); my %codon;my %site;my $type=0;my $s;my $e;my $site;my %type;my %cn;my %gene;my $name;
while(1){
	my $l=<CDS>;
	unless($l){last;}
	chomp $l;
	if(substr($l,0,1) eq ">"){
		if($l=~/complement/){
			$type=1;
			my @a=split"location=complement",$l;
			my @b=split"]",$a[1];
			my @c=split",",$b[0];
			$s=substr($c[0],1);$e=substr($c[1],0,length($c[1])-1);print "$s,$e,$type\n";
		}else{
			$type=0;
			my @a=split"location=",$l;
			my @b=split"]",$a[1];
			my @c=split",",$b[0];###手工将.. 替换成,
			$s=$c[0];$e=$c[1];#print "$s,$e,$type\n";
		}
		$site=0;$name=$l;
	}else{
		unless(substr($l,length($l)-1,1)=~/[0-9a-zA-Z]/){$l=substr($l,0,length($l)-1);}
		my $len=length($l);
		foreach  (1..$len) {
			$site++; my $my;
			if($type==0){
				$my=$s+$site-1;$gene{$my}=$name;#print "$site,$my\n";
				$site{$my}=substr($l,$_-1,1);$type{$my}=$type;
				my $cn=$site%3;$cn{$my}=$cn;
				if($cn==0){
					$codon{$my}=$site{$my-2}.$site{$my-1}.$site{$my}."_3";
					$codon{$my-1}=$site{$my-2}.$site{$my-1}.$site{$my}."_2";
					$codon{$my-2}=$site{$my-2}.$site{$my-1}.$site{$my}."_1";
				}
			}elsif($type==1){
				$my=$e-$site+1;$gene{$my}=$name;#print "$site,$my\n";
				$site{$my}=substr($l,$_-1,1);$type{$my}=$type;
				my $cn=$site % 3;$cn{$my}=$cn;
				if($cn==0){
					$codon{$my}=$site{$my+2}.$site{$my+1}.$site{$my}."_3";print "$site,$my,$codon{$my}\n";
					$codon{$my+1}=$site{$my+2}.$site{$my+1}.$site{$my}."_2";
					$codon{$my+2}=$site{$my+2}.$site{$my+1}.$site{$my}."_1";
				}
			}
		}
	}
}
close CDS;
my @site=sort {$a<=>$b}keys %site;
foreach my $site (@site) {
#	print "$site,$site{$site},$codon{$site},$gene{$site}\n";
}

my(%genetic_code) = (

            'TCA' => 'S',    # Serine
            'TCC' => 'S',    # Serine
            'TCG' => 'S',    # Serine
            'TCT' => 'S',    # Serine
            'TTC' => 'F',    # Phenylalanine
            'TTT' => 'F',    # Phenylalanine
            'TTA' => 'L',    # Leucine
            'TTG' => 'L',    # Leucine
            'TAC' => 'Y',    # Tyrosine
            'TAT' => 'Y',    # Tyrosine
            'TAA' => '',    # Stop
            'TAG' => '',    # Stop
            'TGC' => 'C',    # Cysteine
            'TGT' => 'C',    # Cysteine
            'TGA' => '',    # Stop
            'TGG' => 'W',    # Tryptophan
            'CTA' => 'L',    # Leucine
            'CTC' => 'L',    # Leucine
            'CTG' => 'L',    # Leucine
            'CTT' => 'L',    # Leucine
            'CCA' => 'P',    # Proline
            'CCC' => 'P',    # Proline
            'CCG' => 'P',    # Proline
            'CCT' => 'P',    # Proline
            'CAC' => 'H',    # Histidine
            'CAT' => 'H',    # Histidine
            'CAA' => 'Q',    # Glutamine
            'CAG' => 'Q',    # Glutamine
            'CGA' => 'R',    # Arginine
            'CGC' => 'R',    # Arginine
            'CGG' => 'R',    # Arginine
            'CGT' => 'R',    # Arginine
            'ATA' => 'I',    # Isoleucine
            'ATC' => 'I',    # Isoleucine
            'ATT' => 'I',    # Isoleucine
            'ATG' => 'M',    # Methionine
            'ACA' => 'T',    # Threonine
            'ACC' => 'T',    # Threonine
            'ACG' => 'T',    # Threonine
            'ACT' => 'T',    # Threonine
            'AAC' => 'N',    # Asparagine
            'AAT' => 'N',    # Asparagine
            'AAA' => 'K',    # Lysine
            'AAG' => 'K',    # Lysine
            'AGC' => 'S',    # Serine
            'AGT' => 'S',    # Serine
            'AGA' => 'R',    # Arginine
            'AGG' => 'R',    # Arginine
            'GTA' => 'V',    # Valine
            'GTC' => 'V',    # Valine
            'GTG' => 'V',    # Valine
            'GTT' => 'V',    # Valine
            'GCA' => 'A',    # Alanine
            'GCC' => 'A',    # Alanine
            'GCG' => 'A',    # Alanine
            'GCT' => 'A',    # Alanine
            'GAC' => 'D',    # Aspartic Acid
            'GAT' => 'D',    # Aspartic Acid
            'GAA' => 'E',    # Glutamic Acid
            'GAG' => 'E',    # Glutamic Acid
            'GGA' => 'G',    # Glycine
            'GGC' => 'G',    # Glycine
            'GGG' => 'G',    # Glycine
            'GGT' => 'G',    # Glycine
            );
###
my @file=glob "$file/*/*.snps";
open(O,">$out");
	print O "FILE,Total SNP,SNP in coding,SNP in non-coding,Synonymous,Non-Synonymous\n";
foreach my $f (@file) {
	print "$f\n";
	open(F,$f);
	open(OUT,">$f.circos");
	my $snp;my $csnp;my $ncsnp;my $sy;my $nsy;my %new;
	while(1){
		my $line=<F>;
		unless($line){last;}
		chomp $line;
		my @a=split"\t",$line;
		unless($line=~/NC_018042.1/){next;}
		if($a[1] eq $a[2]){next;}my $tmp=$a[0]+1;
		#print "$line\n";
		unless($a[1]=~/A|T|G|C/ && $a[2]=~/A|T|G|C/){next;}
		$snp++;$new{$a[0]}=$a[2];
		if(exists $codon{$a[0]}){
		my $old_codon=$codon{$a[0]};
		my $n=substr($old_codon,length($old_codon)-1,1);my $new_codon;
		if($n==1){
			$new_codon=$a[2].substr($old_codon,1,2);
		}elsif($n==2){
			if(exists $new{$a[0]-1}){
				$new_codon=$new{$a[0]-1}.$a[2].substr($old_codon,2,1);
			}else{
				$new_codon=substr($old_codon,0,1).$a[2].substr($old_codon,2,1);
			}
		}elsif($n==3){
			my $n1;my $n2;
			if(exists $new{$a[0]-2}){$n1=$new{$a[0]-2};}else{$n1=substr($old_codon,0,1);}
			if(exists $new{$a[0]-1}){$n2=$new{$a[0]-1};}else{$n2=substr($old_codon,1,1);}
			$new_codon=$n1.$n2.$a[2];
		}
		my $old=$genetic_code{substr($old_codon,0,3)};
		my $new=$genetic_code{$new_codon};
		
		print "$a[0],$old_codon,$old,$new_codon,$new\n";
			$csnp++;
			if($old eq $new){$sy++;
			print OUT "Mitochondrion\t$a[0]\t$tmp\tfill_color=purple\n";#同义突变purple
			}else{$nsy++;
			print OUT "Mitochondrion\t$a[0]\t$tmp\tfill_color=blue\n";#非同义突变blue
			}
		}else{
			$ncsnp++;
			print OUT "Mitochondrion\t$a[0]\t$tmp\tfill_color=vlpblack\n";#Non-coding SNP black
		}
	}
	close F;close OUT;
	print "$f,$snp,$csnp,$ncsnp,$sy,$nsy\n";

	print O "$f,$snp,$csnp,$ncsnp,$sy,$nsy\n";
}