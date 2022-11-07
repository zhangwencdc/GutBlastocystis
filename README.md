# GutBlastocystis

## Function
This pipeline is working for identification of Blastocystis sp. subtype 3 in metagenomic samples.
The script is writen by Perl.

## Usage
perl Blastocystis_finder.pl -fq1 <fq1> 
  or
perl Blastocystis_finder.pl -fq1 <fq1>  -fq2 <fq2>

### Input
  fastq format file (SE file or PE files)
### Output
  (1) sam align files
  (2) snp list to ref seq
  (3) Report files (Covereage,Depth)
  
