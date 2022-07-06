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
		#print "/home/zhangwen/bin/bowtie2-2.4.4-linux-x86_64/bowtie2 -1 $fq1 -2 $fq2 -x $ref -S $a[0].target.sam --no-unal\n";
system "/home/zhangwen/bin/bowtie2-2.4.4-linux-x86_64/bowtie2 -1 $fq1 -2 $fq2 -x $ref -S $a[0].target.sam --no-unal\n";
system "samtools view -bS $a[0].target.sam > $a[0].target.bam\n";
system "samtools sort $a[0].target.bam -o $a[0].target.sort.bam\n";
system "samtools index $a[0].target.sort.bam\n";
system "bcftools mpileup -Ou -f  $ref $a[0].target.sort.bam | bcftools call -m -Ov -o  $a[0].target.vcf\n";
	
		###
		open(F,"$a[0].target.sam");
		#print "$fq1,$fq2,$o\n";
		#print "$a[0].virus.sam\n";
		my %read;
		while(1){
			my $line=<F>;
			unless($line){last;}
			chomp $line;
			if(substr($line,0,1) eq "@"){next;}
			my @a=split" ",$line;
			unless($a[2]=~/[0-9a-zA-Z]/){next;}
			my $len=length($a[9]);
			
			my $tmp;
			foreach my $a (@a) {
				if($a=~/^MD:Z:/){$tmp=$a;}
			}
			my $l=length($tmp);
			my $num="";my $match=0;my $total=$len;
			foreach  (0..($l-1)) {
				my $site=substr($tmp,$_,1);
				if($site=~/[0-9]/){
					$num=$num.$site;
				}else{
					$match+=$num;

					$num="";			
				}
			}
			$match+=$num;
			
			unless($len>0){next;}
			unless($match>0){next;}
			my $maper=$match/$len*100;
			if($match<$match_cutoff && $maper<$cutoff){next;}  ###
			my @b=split"/",$a[0];
			$read{$b[0]}++;


		}
		close F;


		my @key=keys %read; my %virus;
		my %result; my $num=0;
		open(OL,">$a[0].matchread");
		foreach my $key (@key) {
			if($read{$key}<1){next;} #没有与virus比对
			
			$result{$key}++;
			print OL "$key\n";
			$num++;
		}
		close OL;
		print "Candidate paired reads uniq map to Target is $num\n";

		   open(O1,">$a[0].fq1.filter.fq");
		open(F1,$fq1);my $abs;
		while(1){
				my $line=<F1>;
				unless($line){last;}
				chomp $line;
				my $seq=<F1>;
				chomp $seq;
				my $a=<F1>;
				chomp $a;
				my $b=<F1>;
				chomp $b;
				$line=substr($line,1);
				my @name=split" ",$line;
				my @b=split"/",$name[0];
				my $name=$b[0];

			  # print "$name,$read{$name}\n";
				$abs++;
				if(exists $read{$name}){
					
						print O1 "@";
						print O1 "$name\n$seq\n$a\n$b\n";
				}

		}

		close F1;close O1;

		if(defined $fq2){
			open(O2,">$a[0].fq2.filter.fq");
		open(F2,$fq2);
		while(1){
				my $line=<F2>;
				unless($line){last;}
				chomp $line;
				my $seq=<F2>;
				chomp $seq;
				my $a=<F2>;
				chomp $a;
				my $b=<F2>;
				chomp $b;
				$line=substr($line,1);
				my @name=split" ",$line;
				my @b=split"/",$name[0];
				my $name=$b[0];

		#       print "$name\n";
				if(exists $result{$name}){
					
						print O2 "@";
						print O2 "$name\n$seq\n$a\n$b\n";
				}

		}
		close F2;
		}
		close O2;
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
				
			}elsif(length($b[3])==1 && length($b[4])==1 && $b[5]<50 && $dp>=10){
				if($b[3] ne $b[4]){print OH "$a[0] $line\n";$he++;}
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