#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
my $outdir||=getcwd;
my $threshold||=0.00005;
my (@factor,$qiime,@input,$format);
GetOptions(
	'help|?'=>sub{usage{}},
	'i=s{1,}'=>\@input,
	'f:s'=>\$format,
	'factor=s{1,}'=>\@factor,
	'o:s'=>\$outdir,
	'p:s'=>\$threshold,
	'qiime=s'=>\$qiime
);

if (!@input || !$format) {
	&usage();
}

open(OUT,">$qiime");
for (my $k=0;$k<=$#input;$k++)
{
	my $file = -e "$outdir/temp.fasta";
	if ($file)
	{
		system "rm $outdir/temp.fasta";
	}
	my $i=0;
	if ($format =~ /fastq/ || $format =~/fq/)
	{
		system "/share/work/biosoft/fastx_toolkit/v0.0.13/fastq_to_fasta -i $input[$k] -o  $outdir/temp.fasta";	
	}
	if ($format =~ /fasta/ || $format =~/fa/)
	{
		system "cp $input[$k] $outdir/temp.fasta";
	}
	open(IN,"temp.fasta");
	while(<IN>)
	{
		chomp;
		if ($_=~/\>/)
		{
			$i++;
			print OUT ">",$factor[$k],"_",$i,"\n";
		}
		else
		{
			print OUT $_,"\n";
		}
	}
}
my $file = -e "$outdir/temp.fasta";
if ($file)
{
	system "rm $outdir/temp.fasta";
}


my $seq_num = (split /\s+/, `wc -l $qiime`)[0];
my $filter_num = int(($seq_num/2)*$threshold)+1;


sub usage{
	print qq{
	
This script deal with fastq file, output could be used by qiime.

Optionts:

 -i			input file
 
 -f			the format of input(fq(fastq) or fa(fasta))

 -factor	the factor of sample

 -qiime		output file
 
 -help		print a brief help message and exits

for example:
perl $0 -i BC1.fastq BC2.fastq -f fq -factor BC1 BC2 -qiime qiime.fasta

or

perl $0 -i BC1.fasta BC2.fasta -f fa -factor BC1 BC2 -qiime qiime.fasta

Email:fanyucai1\@126.com
vesion1.0
2015.9.25		
	};
	exit;
}