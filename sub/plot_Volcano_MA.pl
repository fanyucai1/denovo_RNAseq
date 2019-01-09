#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd;
use Getopt::Long;

my $R="/share/work/biosoft/R/R-v3.0.0/bin/";
my($input,$outdir,$vs);
my $logFC||=1;
my $FDR||=0.001;
my $png_height||=3000;
my $png_width||=3000;
GetOptions(
    "i:s"=>\$input,
    "vs:s"=>\$vs,
    "o:s"=>\$outdir,
    "logFC:s"=>\$logFC,
    "height:s"=>\$png_height,
    "width:s"=>\$png_width,
           );
sub usage{
    print qq{
This script plot Volcano and MA plot.
usage:
perl $0 -i matrix.txt.cond_A_vs_cond_B.edgeR.DE_results -vs A_vs_B -o /path/to/directory/
options:
-i          the input from edgeR output,the file consist 5 cloum:"ID   logFC   logCPM  PValue  FDR"(force)
-vs         the prefix of output(force)
-o          output directory(force)
-height     the picture of height,default:3000
-width      the picture of width,default:3000
logFC       default:1
FDR         default:0.001
};
    exit;
}

if (!$input||!$outdir) {
    &usage();
}



open(IN, $input) ;
open(OUT,">$outdir/$vs.xls");
while (<IN>)
{
    chomp;
    if ($_=~/logFC/)
    {
      print OUT "ID\tlogFC\tlogCPM\tPValue\tFDR\tregulate\n"; 
    }
    else
    {
        my @array=split;
        if ($array[1]>$logFC)
        {
            print OUT "$_\tup\n";
        }
        elsif ($array[1]<(-$logFC))
        {
            print OUT "$_\tdown\n";
        }
        else
        {
            print OUT "$_\tnoraml\n";
        }   
    }
}

system "echo '
#!$R/Rscript
library(ggplot2)
a<-read.table(\"$outdir/$vs.xls\",header=T)
#Volcano plot 
png(filename=\"$outdir/$vs\_Volcano.png\", height = $png_height, width = $png_width, res = $png_height/6, units = \"px\")
p<-ggplot(data=a,aes(x=logFC,y=-log10(FDR),colour=factor(regulate)))+geom_point(size=1)+ scale_colour_manual(values=c(\"green\",\"black\",\"red\"))
p<-p+geom_vline(xintercept=c(-$logFC,$logFC), linetype=\"longdash\", size=0.2)
p<-p + geom_hline(yintercept=c(-log10($FDR)), linetype=\"longdash\", size=0.2)
p<-p + labs(list(title=\"$vs\_Volcano\", x=\"log2(FC)\"))
print(p)
dev.off()
#MA plot
png(filename=\"$outdir/$vs\_MA.png\", height = $png_height, width = $png_width, res = $png_height/6, units = \"px\")
p<-ggplot(data=a,aes(x=logCPM,y=logFC,colour=factor(regulate)))+geom_point(size=1)+ scale_colour_manual(values=c(\"green\",\"black\",\"red\"))
p<-p + labs(list(title=\"$vs\_MA\", x=\"logCPM\"))
print(p)
dev.off()

'>$outdir/draw_MA_Volcano.R";


system "$R/Rscript $outdir/draw_MA_Volcano.R";






