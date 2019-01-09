#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin);

my $fasta;
my $bins ||=100;
my $outdir||=getcwd;
my $min||=200;
my $max||=3000;
my $pdf_width||=15;
my $pdf_height||=7;
my $png_width||=1200;
my $png_height||=700;
my $R="/share/work/biosoft/R/R-v3.0.0/bin/Rscript";
my $xlab||="Length(bp)";
my $ylab||="Number of unigene";

GetOptions(
    "max:s"=>\$max,
    "min:s"=>\$min,
    "bins:s"=>\$bins,
    "fasta:s"=>\$fasta,
    "o:s"=>\$outdir,
    "pngw:s"=>\$png_width,
    "pngh:s"=>\$png_height,
    "pdfw:s"=>\$pdf_width,
    "pdfh:s"=>\$pdf_height,
    "xlab:s"=>\$xlab,
    "ylab:s"=>\$ylab,
           );

sub usage{
    print qq{
This script could plot fasta length distribution.
usage:
perl $0 -fasta file.fa -o /path/to/directory/
options:
-fasta      input fasta file
-o          the output directory
-max        the max length of plot(default:3000)
-min        the mix length of plot(default:200)
-bins       the interval(default:100)
-pdfh       height of pdf(default:7)
-pdfw       width of pdf(default:12)
-pngw       width of png(default:900)
-pngh       height of png(default:700)
-xlab       xlab
-ylab       ylab
    };
    exit;
}

if (!$fasta) {
    &usage();
}

open(IN, $fasta);
open(OUT,">$outdir/length.txt");    
my $seqname;
my %hash;
while (<IN>)
{
    chomp;
    if ($_=~/\>/)
    {
        $seqname=substr($_,1);
    }
    else
    {
        $hash{$seqname}.=$_;
    }  
}
print OUT "length\tcount\n";
my %hash2;
foreach my $key(keys %hash)
{
    if (int(length($hash{$key})/$bins)*$bins<$max)
    {
        $hash2{int(length($hash{$key})/$bins)*$bins}++;
    }
    else
    {
        $hash2{$max}++;
    }
}
for (my $k=$min;$k<=$max;$k=$k+$bins)
{
    if ($k==$max)
    {
        print OUT ">$k","\t",$hash2{$k},"\n";
    }
    else
    {
        my $num1=$k+1;
        my $num2=$k+100;
        print OUT $num1,"-",$num2,"\t",$hash2{$k},"\n";
    }
    
}


system "echo '
#!$R
a<-read.table(file=\"$outdir/length.txt\",header = T)
b<- data.frame(a)
library(ggplot2)
png(\"$outdir/Unigene_length.png\", width=$png_width, height=$png_height)
ggplot(data=b,aes(factor(length,order=T,levels = length),count,label=count))+
  geom_bar(stat=\"identity\",fill=\"#DD8888\",width=0.5, position = \"dodge\")+
  xlab(\"$xlab\")+ylab(\"$ylab\")+
  geom_text(aes(count = count+5), position = position_dodge(0.9), vjust = 0)+
  theme(panel.grid =element_blank(),axis.title.x = element_text(face=\"bold\", size=15),axis.text.x  = element_text(angle=45, vjust=0.5, size=10),axis.title.y = element_text(face=\"bold\", size=15))
dev.off()
pdf(\"$outdir/Unigene_length.pdf\", width=$pdf_width, height=$pdf_height)
ggplot(data=b,aes(factor(length,order=T,levels = length),count,label=count))+
  geom_bar(stat=\"identity\",fill=\"#DD8888\",width=0.5, position = \"dodge\")+
  xlab(\"$xlab\")+ylab(\"$ylab\")+
  geom_text(aes(count = count+5), position = position_dodge(0.9), vjust = 0)+
  theme(panel.grid =element_blank(),axis.title.x = element_text(face=\"bold\", size=15),axis.text.x  = element_text(angle=45, vjust=0.5, size=10),axis.title.y = element_text(face=\"bold\", size=15))
dev.off()
'>$outdir/length_hist.R";

system "$R $outdir/length_hist.R";
