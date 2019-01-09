#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
my $R="/share/work/biosoft/R/R-v3.2.3/bin/";
my ($FPKM,$outdir);
GetOptions(
    "i:s"=>\$FPKM,
    "o:s"=>\$outdir, 
           );

sub usage{
    print qq{
This script plot the boxplot and density of FPKM.The FPKM file is from Trinity(named:####.TMM_normalized.FPKM).
usage:
perl $0 -i matrix.txt.TMM_normalized.FPKM -o /path/to/outdirectoy
    };
    exit;
}

if (!$FPKM  || !$outdir) {
    &usage();
}


my $num=`awk '{print NF}' $FPKM |tail -n1`;
chomp($num);
open(SH,">$outdir/temp.sh");
system "echo 'FPKM\tsample'>$outdir/FPKM_boxplot_density_input.txt";
for(my $i=2;$i<=$num;$i++)
{
    my $j=$i-1;
    my $name=`awk '{print \$$j}' $FPKM |head -n1`;
    chomp($name);
    print SH "sed 1d $FPKM|awk \'{print \$$i\"\\t$name\"}\'>>$outdir/FPKM_boxplot_density_input.txt\n";#code
}
system "sh $outdir/temp.sh\n";
system "echo '
#!$R/Rscript
library(ggplot2)
a=read.table(\"$outdir/FPKM_boxplot_density_input.txt\",header = T)
png(filename=\"$outdir/FPKM_density.png\",height=800,width=1000);
ggplot(a,aes(log10(FPKM),colour=sample))+geom_density()+theme(axis.title=element_text(size=20),axis.text=element_text(size=20))
dev.off()
png(filename=\"$outdir/FPKM_boxplot.png\",height=800,width=1500);
ggplot(a,aes(sample,log10(FPKM),fill=sample))+geom_boxplot()+theme(axis.title=element_text(size=20),axis.text=element_text(size=20))
dev.off()
pdf(\"$outdir/FPKM_density.pdf\");
ggplot(a,aes(log10(FPKM),colour=sample))+geom_density()
dev.off()
pdf(\"$outdir/FPKM_boxplot.pdf\");
ggplot(a,aes(sample,log10(FPKM),fill=sample))+geom_boxplot()
dev.off()
'>$outdir/FPKM_boxplot_density.Rscript";

system "$R/Rscript $outdir/FPKM_boxplot_density.Rscript\n";