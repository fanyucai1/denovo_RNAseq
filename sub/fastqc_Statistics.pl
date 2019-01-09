#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd;
use Getopt::Long;
my ($pe1,$pe2,$outdir,$sampleID,$prefix);
GetOptions(
     "a:s"=>\$pe1,       
     "b:s"=>\$pe2,
     "o:s"=>\$outdir,
     "id:s"=>\$sampleID,
     "p:s"=>\$prefix,
           );
my $fastq_qc_stat="/share/work/biosoft/fastq_qc_stat/fastq_qc_stat";
sub usage{
    print qq{
This script could satistics the PE_fastq file.
usage:
perl $0 -a sample1_1.fq -b sample1_2.fq /path/to/diretory -p prefix
options:
-a        input fastq file(force)
-b        input fastq file
-o        the output directory(force)
-p        the prefix of output
    };
    exit;
}
if (!$pe1 || !$outdir||!$sampleID) {
    &usage();
}
my @id=split(/\,/,$sampleID);
my @left=split(/\,/,$pe1);
if($#id>0)
{
   open(OUT,">$outdir/statistics_fq.table");
}
else
{
   open(OUT,">$outdir/$prefix.statistics_fq.table");
}
print OUT "Samples\tRead_Number\tBase_Number\tGC_Content(%)\t%>Q20\t%>Q30\n";
if (defined $pe2)
{
    my @right=split(/\,/,$pe2);
    for(my $k=0;$k<=$#left;$k++)
    {
        my $num=0;
        `export LD_LIBRARY_PATH=/share/work/biosoft/R/R-vv3.2.3/lib64/R/lib/:/share/work/biosoft/lib/:\$LD_LIBRARY_PATH`;
        system `$fastq_qc_stat -a $left[$k] -b $right[$k] -f $outdir/$id[$k]`;
        open(IN, "$outdir/$id[$k].stat") ;
        while (<IN>)
        {
            chomp;
            my @array=split;
            if ($_!~/\#/)
            {
                ++$num;
                if ($num==3)
                {
                    print OUT "$id[$k]\t$array[1]\t$array[2]\t$array[3]\t$array[5]\t$array[7]\n";
                }
            } 
        }
    }   
}
close;
