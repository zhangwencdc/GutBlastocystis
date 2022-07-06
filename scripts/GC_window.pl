#!/usr/bin/perl

my $file=$ARGV[0];#input
my $window=$ARGV[1];#100
my $step=$ARGV[2];#10

open(FILE,"$file")||die;
open(OUT,">GC_window.txt");

my $value=0; my $seq;
while(1){
        my $line=<FILE>;
        unless($line){
                        last;
        }
        chomp $line;
		if($line=~/>/){next;}
		unless(substr($line,length($line)-1,1)=~/[a-zA-Z]/){$line=substr($line,0,length($line)-1);}
		$seq.=$line;
      }
	  my $n=int(length($seq)/$step);
	  foreach  (0..$n) {
		  my $start=$_*$step;
		  my $end=$start+$window-1;
		  my $s=substr($seq,$start,$window);
		  my $gc=&gc($s);
		  print OUT "Mitochondrion\t$start\t$end\t$gc\n";
	  }
#################
sub gc{
my $seq=shift @_;
my $gc=$seq=~tr/(G|C|g|c)/(G|C|g|c)/;
my $l=length($seq);
my $g=$gc/$l;

return $g;
}
