#!/usr/bin/perl

=head1 Name:

=head1 Description

Identification for Blastocystis 

=head1 Version
        Version: 1.0,  Date: 2022-09-26

=head1 Usage

  perl  Blastocystisis.pl [options]

  --keyname     Tag; 默认为Input
  --fq1         Fq1
  --fq2                 Fq2
  --outdir      输出路径 默认为当前路径
   --type
  --help
=head1 Example
perl Blastocystisis.pl -Key test -fq1 fq1.fastq -fq2 fq2.fastq

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin);
use File::Basename qw(basename dirname);
use Data::Dumper;
use File::Path;
use Cwd;
my $path = getcwd;

##
my ($Keyname,$Ref,$fq1,$fq2,$database,$type,$fa,$cutoff,$Outdir,$verbose,$db,$match_cutoff,$cutoff);
my ($Verbose,$Help);
#my $Time_Start = sub_format_datetime(localtime(time())); #......
#my $Data_Vision = substr($Time_Start,0,10);

GetOptions(
      
        "fq1|1:s"=>\$fq1,
                "fq2|2:s"=>\$fq2,
        "cutoff:s"=>\$cutoff,
        "outdir:s"=>\$Outdir,
        "verbose"=>\$Verbose,
        "Key|K|Keyname:s" =>\$Keyname,
		"m|ML:s" => \$match_cutoff,
		"c|C:s" => \$cutoff,
        "help"=>\$Help
);
my $Time_Start= sub_format_datetime(localtime(time()));
if(!defined $Outdir){$Outdir="./";}
if ( !defined $match_cutoff ){ $match_cutoff=70;}
if ( !defined $cutoff ){ $cutoff=80; }
if(!defined $Keyname){$Keyname="Input";}


my $bowtie2build="bowtie2-build";
my $bowite="bowtie2";
my $target=$Bin."/Blastocystis-sub3_Genome.fasta";
my $o=$Outdir."/".$Keyname."_G";
unless(-e "$target.1.bt2"){system "$bowtie2build $target $target\n";}

open(FILE,"$target");
my %gene;my %genelen;my $len;my $name;
while(1){
        my $line=<FILE>;
        unless($line){$genelen{$name}=$len;last;}
        chomp $line;
        unless(substr($line,length($line)-1,1)=~/[0-9a-zA-Z]/){$line=substr($line,0,length($line)-1);}
        unless($line=~/>/){$len+=length($line);next;}
        $genelen{$name}=$len;
   #     print "$name,$len\n";
        my @a=split" ",$line;
        $name=substr($a[0],1);$gene{$name}=$line;


        $len=0;
}
close FILE;


print "Step1: Target Genome align\n";
if(defined $fq2){system "$bowite -1 $fq1 -2 $fq2 -x $target -S $o.genome.sam --no-unal\n";}else{system "$bowite -U $fq1  -x $target -S $o.genome.sam --no-unal\n";}
system "samtools view -bS $o.genome.sam > $o.genome.bam\n";
system "samtools sort $o.genome.bam -o $o.genome.sort.bam\n";
system "samtools index $o.genome.sort.bam\n";
system "bcftools mpileup -Ou -f  $target $o.genome.sort.bam | bcftools call -mv -Ov -o  $o.genome.vcf\n";

###
open(FILE,"$o.genome.sam");
#print "$fq1,$fq2,$o\n";
#print "$o.virus.sam\n";
my %read;
while(1){
        my $line=<FILE>;
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
close FILE;


my @key=keys %read; my %virus;
my %result; my $num=0;
open(OL,">$o.matchread");
foreach my $key (@key) {
        if($read{$key}<1){next;} 

        $result{$key}++;
        print OL "$key\n";
        $num++;
}
close OL;
print "Candidate paired reads uniq map to Target is $num\n";

   open(OUT,">$o.fq1.filter.fasta"); open(O1,">$o.fq1.filter.fq");
open(FILE,$fq1);my $abs;
while(1){
        my $line=<FILE>;
        unless($line){last;}
        chomp $line;
        my $seq=<FILE>;
        chomp $seq;
        my $a=<FILE>;
        chomp $a;
        my $b=<FILE>;
        chomp $b;
        $line=substr($line,1);
        my @name=split" ",$line;
        my @b=split"/",$name[0];
        my $name=$b[0];

      # print "$name,$read{$name}\n";
                $abs++;
        if(exists $read{$name}){
               print OUT ">";
                print OUT "$name fq1\n$seq\n";
                                print O1 "@";
                                print O1 "$name\n$seq\n$a\n$b\n";
        }

}

close FILE;close OUT;

if(defined $fq2){
        open(OUT,">$o.fq2.filter.fasta");open(O2,">$o.fq2.filter.fq");
open(FILE,$fq2);
while(1){
        my $line=<FILE>;
        unless($line){last;}
        chomp $line;
        my $seq=<FILE>;
        chomp $seq;
        my $a=<FILE>;
        chomp $a;
        my $b=<FILE>;
        chomp $b;
        $line=substr($line,1);
        my @name=split" ",$line;
        my @b=split"/",$name[0];
        my $name=$b[0];

#       print "$name\n";
        if(exists $result{$name}){
                        print OUT ">";
                print OUT "$name fq2\n$seq\n";
                                print O2 "@";
                                print O2 "$name\n$seq\n$a\n$b\n";
        }

}
close FILE;
}
close OUT;
print "Step2: Target Genome filter align\n";
if(defined $fq2){system "$bowite -1 $o.fq1.filter.fasta -2 $o.fq2.filter.fasta -x $target -S $o.genome.filter.sam -f --no-unal\n";}else{system "$bowite -U $o.fq1.filter.fasta  -x $target -S $o.genome.filter.sam -f --no-unal\n";}
open(FILE,"$o.genome.filter.sam");
#print "$fq1,$fq2,$o\n";
#print "$o.virus.sam\n";
my %read;my %start;my %end;my %align;
while(1){
        my $line=<FILE>;
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
        $read{$a[2]}++;
        if($a[1] == 16){ 
                        foreach  (($a[3]-length($a[9])+1)..$a[3]) {
                                $align{$a[2]}{$_}++;
                }

        }else{

                foreach  ($a[3]..($a[3]+length($a[9])-1)) {
                                $align{$a[2]}{$_}++;
                }
        }


}
close FILE;


open(OUT,">$o.genome.report");
print OUT "Gene ID,Target Length,Match Read Num,Align Length,Coverage%,Depth,Gene Annotation,\n";
my @key=sort {$read{$b}<=>$read{$a}} keys %read; my $sum;
foreach my $key (@key) {
        print OUT "$key,$genelen{$key},$read{$key},";my $align;my $depth;
        foreach  (1..$genelen{$key}) {
                if($align{$key}{$_}>=1){$align++;$depth+=$align{$key}{$_};}
        }
        if($align>0){$depth=$depth/$align;}
        my $coverage=$align/$genelen{$key}*100;
        $coverage=sprintf "%.2f",$coverage;
        print OUT "$align,$coverage %,$depth,$gene{$key}\n";
		if($align>=1000){$sum+=$align;}
}
close OUT;

my $target=$Bin."/Blastocystis_sub3_Mitochondria.fasta";
my $o=$Outdir."/".$Keyname."_M";
unless(-e "$target.1.bt2"){system "$bowtie2build $target $target\n";}

open(FILE,"$target");
my %gene;my %genelen;my $len;my $name;
while(1){
        my $line=<FILE>;
        unless($line){$genelen{$name}=$len;print "$name,$len\n";last;}
        chomp $line;
        unless(substr($line,length($line)-1,1)=~/[0-9a-zA-Z]/){$line=substr($line,0,length($line)-1);}
        unless($line=~/>/){$len+=length($line);next;}
        $genelen{$name}=$len;
   #     print "$name,$len\n";
        my @a=split" ",$line;
        $name=substr($a[0],1);$gene{$name}=$line;


        $len=0;
}
close FILE;


print "Step1: Target Mitochondria align\n";
if(defined $fq2){system "$bowite -1 $fq1 -2 $fq2 -x $target -S $o.Mitochondria.sam --no-unal\n";}else{system "$bowite -U $fq1  -x $target -S $o.Mitochondria.sam --no-unal\n";}
system "samtools view -bS $o.Mitochondria.sam > $o.Mitochondria.bam\n";
system "samtools sort $o.Mitochondria.bam -o $o.Mitochondria.sort.bam\n";
system "samtools index $o.Mitochondria.sort.bam\n";
system "bcftools mpileup -Ou -f  $target $o.Mitochondria.sort.bam | bcftools call -mv -Ov -o  $o.Mitochondria.vcf\n";

###
open(FILE,"$o.Mitochondria.sam");
#print "$fq1,$fq2,$o\n";
#print "$o.virus.sam\n";
my %read;
while(1){
        my $line=<FILE>;
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
close FILE;


my @key=keys %read; my %virus;
my %result; my $num=0;
open(OL,">$o.matchread");
foreach my $key (@key) {
        if($read{$key}<1){next;} 

        $result{$key}++;
        print OL "$key\n";
        $num++;
}
close OL;
print "Candidate paired reads uniq map to Target is $num\n";

   open(OUT,">$o.fq1.filter.fasta"); open(O1,">$o.fq1.filter.fq");
open(FILE,$fq1);my $abs;
while(1){
        my $line=<FILE>;
        unless($line){last;}
        chomp $line;
        my $seq=<FILE>;
        chomp $seq;
        my $a=<FILE>;
        chomp $a;
        my $b=<FILE>;
        chomp $b;
        $line=substr($line,1);
        my @name=split" ",$line;
        my @b=split"/",$name[0];
        my $name=$b[0];

      # print "$name,$read{$name}\n";
                $abs++;
        if(exists $read{$name}){
               print OUT ">";
                print OUT "$name fq1\n$seq\n";
                                print O1 "@";
                                print O1 "$name\n$seq\n$a\n$b\n";
        }

}

close FILE;close OUT;

if(defined $fq2){
        open(OUT,">$o.fq2.filter.fasta");open(O2,">$o.fq2.filter.fq");
open(FILE,$fq2);
while(1){
        my $line=<FILE>;
        unless($line){last;}
        chomp $line;
        my $seq=<FILE>;
        chomp $seq;
        my $a=<FILE>;
        chomp $a;
        my $b=<FILE>;
        chomp $b;
        $line=substr($line,1);
        my @name=split" ",$line;
        my @b=split"/",$name[0];
        my $name=$b[0];

#       print "$name\n";
        if(exists $result{$name}){
                        print OUT ">";
                print OUT "$name fq2\n$seq\n";
                                print O2 "@";
                                print O2 "$name\n$seq\n$a\n$b\n";
        }

}
close FILE;
}
close OUT;
print "Step2: Target Mitochondria filter align\n";
if(defined $fq2){system "$bowite -1 $o.fq1.filter.fasta -2 $o.fq2.filter.fasta -x $target -S $o.Mitochondria.filter.sam -f --no-unal\n";}else{system "$bowite -U $o.fq1.filter.fasta  -x $target -S $o.Mitochondria.filter.sam -f --no-unal\n";}
open(FILE,"$o.Mitochondria.filter.sam");
#print "$fq1,$fq2,$o\n";
#print "$o.virus.sam\n";
my %read;my %start;my %end;my %align;
while(1){
        my $line=<FILE>;
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
        $read{$a[2]}++;
        if($a[1] == 16){ 
                        foreach  (($a[3]-length($a[9])+1)..$a[3]) {
                                $align{$a[2]}{$_}++;
                }

        }else{

                foreach  ($a[3]..($a[3]+length($a[9])-1)) {
                                $align{$a[2]}{$_}++;
                }
        }


}
close FILE;


open(OUT,">$o.Mitochondria.report");
print OUT "Gene ID,Target Length,Match Read Num,Align Length,Coverage%,Depth,Gene Annotation,\n";
my @key=sort {$read{$b}<=>$read{$a}} keys %read; my $sum2;
foreach my $key (@key) {
        print OUT "$key,$genelen{$key},$read{$key},";my $align;my $depth;
        foreach  (1..$genelen{$key}) {
                if($align{$key}{$_}>=1){$align++;$depth+=$align{$key}{$_};}
        }
        if($align>0){$depth=$depth/$align;}
        my $coverage=$align/$genelen{$key}*100;
        $coverage=sprintf "%.2f",$coverage;
        print OUT "$align,$coverage %,$depth,$gene{$key}\n";
		if($align>=1000){$sum2+=$align;}
}
close OUT;
print "Genome Align $sum\n Mitochondira Align $sum2\n";
if($sum>=2000000 || $sum2>=0.8*28269){print "$Keyname,Positive\n";}

my $Time_End= sub_format_datetime(localtime(time()));
print "Running from [$Time_Start] to [$Time_End]\n";



#====================================================================================================================
#  +------------------+
#  |   subprogram     |
#  +------------------+



sub sub_format_datetime #.....
{
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
        $wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon, $day, $hour, $min, $sec);
}
#####




