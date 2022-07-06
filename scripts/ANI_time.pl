#!/usr/bin/perl
use strict;
use warnings;

#ANI and time

my $file=$ARGV[0];#CHMP1.ANI.list
my $info=$ARGV[1];#CHMP1.Mito.list

open(FILE,$info);my %time;
while(1){
	my $line=<FILE>;
	unless($line){last;}
	chomp $line;
	my @a=split"\t",$line;
	$time{$a[0]}=$a[1];
}

close FILE;

open(FILE,$file);
open(OUT,">$file.time");
while(1){
	my $line=<FILE>;
	unless($line){last;}
	chomp $line;
	my @a=split"\t",$line;
	my $a=substr($a[0],0,length($a[0])-14);
	my $b=substr($a[1],0,length($a[1])-14);
	my $t1=$time{$a};
	my $t2=$time{$b};
	print OUT "$line\t$a\t$b\t$t1\t$t2\t";
	$t1=substr($t1,3);$t2=substr($t2,3);print "$a,$b,$t1,$t2\n";
	my $dis=$t2-$t1;
	if($dis<0){$dis=0-$dis;}
	print OUT "$dis\n";
}
close FILE;