#!/usr/bin/perl
use strict;
use warnings;
use File::Basename qw(basename dirname);


my @file=glob "/home/zhangwen/project/2022Time/WGS_Analysis/Eukaryote/New_Blastocystis/Blastocystis_Mitochondria/seq/*.fas";

foreach my $f1(@file){
	my $f1n=basename($f1);
	foreach my $f2 (@file) {
		if($f1 eq $f2){next;}
		my $f2n=basename($f2);my $name=$f1n."_".$f2n;
		system "dnadiff $f1 $f2 -p $name\n";
	}
}