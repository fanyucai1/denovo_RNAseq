#!/usr/bin/perl -w
use strict;
use warnings;
use  Cwd;
use Getopt::Long;
my $R="/share/work/biosoft/R/R-v3.0.0/bin";
my($SSR,$outdir);
GetOptions(
    "i:s"=>\$SSR,
    "o:s"=>\$outdir,
        );

sub usage{
    print qq{
This script will plot the SSR result from MISA output.
usage:
perl $0 -i unigene.fasta.statistics -o /path/to/directory
options:
-i              input SSR result named:###.statistics
-o              the output directory
    };
    exit;
}
if (!$SSR || !$outdir) {
    &usage();
}

open(IN,"$SSR");
my %hash;
my $num=0;
my %p;
while (<IN>)
{
    chomp;
    if ($_=~/Frequency of classified repeat types/)
    {
        ++$num;
    }
    if ($num==1)
    {
        if ($_=~/[ACGT]/ && $_=~/\//)
        {
            my @array=split(/\t/,$_);
            for(my $k=1;$k<=$#array;$k++)
            {
                if ($array[$k]=~/[0-9]/)
                {
                    $hash{$array[0]}+=$array[$k];
                }  
            }  
        }    
    }    
}
foreach my $key (keys %hash)
{
    $hash{$key}=$hash{$key}/2;
    if (length($key)==3)
    {
        $p{"Mono nucleotide"}+=$hash{$key};
    }
    if(length($key)==5)
    {
        $p{"Di nucleotide"}+=$hash{$key};
    }
    if(length($key)==7)
    {
        $p{"Tri nucleotide"}+=$hash{$key};
    }
    if(length($key)==9)
    {
        $p{"Tetra nucleotide"}+=$hash{$key};
    }
    if(length($key)==11)
    {
        $p{"Penta nucleotide"}+=$hash{$key};
    }
    if(length($key)==13)
    {
        $p{"Hexa nucleotide"}+=$hash{$key};
    }
}
open(OUT,">$outdir/SSR_plot.txt");
print OUT "type\tnumber\tkind\n";
print OUT "Mono_nucleotide\t$p{'Mono nucleotide'}\t1\n";
foreach my $key (keys %hash)
{
    my $len=length($key);
    if ($len == 3)
    {
        print OUT "$key\t$hash{$key}\t1\n";
    }  
}
print OUT "Di_nucleotide\t$p{'Di nucleotide'}\t2\n";
foreach my $key (keys %hash)
{
     my $len=length($key);
     if ($len== 5)
    {
        print OUT "$key\t$hash{$key}\t2\n";
    }  
}
print OUT "Tri_nucleotide\t$p{'Tri nucleotide'}\t3\n";
foreach my $key (keys %hash)
{
     my $len=length($key);
     if ($len== 7)
    {
        print OUT "$key\t$hash{$key}\t3\n";
    }  
}
print OUT "Tetra_nucleotide\t$p{'Tetra nucleotide'}\t4\n";
foreach my $key (keys %hash)
{
     my $len=length($key);
     if ($len== 9)
    {
        print OUT "$key\t$hash{$key}\t4\n";
    }  
}
print OUT "Penta_nucleotide\t$p{'Penta nucleotide'}\t5\n";
foreach my $key (keys %hash)
{
     my $len=length($key);
    if ($len== 11)
    {
        print OUT "$key\t$hash{$key}\t5\n";
    }  
}
print OUT "Hexa_nucleotide\t$p{'Hexa nucleotide'}\t6\n";
foreach my $key (keys %hash)
{
     my $len=length($key);
    if ($len== 13)
    {
        print OUT "$key\t$hash{$key}\t6\n";
    }  
}

system "echo '
#!$R/Rscript
library(ggplot2)
a<-read.table(\"$outdir/SSR_plot.txt\",header = T)
png(filename = \"$outdir/SSR_Plot.png\",height = 2000,width=5500)
p<-ggplot(a,aes(factor(type,order=T,levels = type),number,label=number))+geom_bar(stat = \"identity\",width = 0.5,aes(fill=as.character(kind)))+scale_fill_brewer(palette=\"Set1\")+theme(legend.position=\"none\")
p+geom_text(aes(number = number + 0.05), position = position_dodge(0.9), vjust = 0,size=9)+labs(list(x=\"Repeat nucleotide types\",y=\"SSR motif numbers\"))+theme(axis.title = element_text(face=\"bold\", colour=\"#990000\", size=25),axis.text  = element_text(angle=30, vjust=0.5, size=20))
dev.off()
pdf(\"$outdir/SSR_Plot.pdf\",width=70,height=10)
p<-ggplot(a,aes(factor(type,order=T,levels = type),number,label=number))+geom_bar(stat = \"identity\",width = 0.5,aes(fill=as.character(kind)))+scale_fill_brewer(palette=\"Set1\")+theme(legend.position=\"none\")
p+geom_text(aes(number = number + 0.05), position = position_dodge(0.8), vjust = 0,size=9)+labs(list(x=\"Repeat nucleotide types\",y=\"SSR motif numbers\"))+theme(axis.title = element_text(face=\"bold\", colour=\"#990000\", size=25),axis.text  = element_text(angle=30, vjust=0.5, size=20))
dev.off()
'>$outdir/SSR.Rscript";

system "$R/Rscript $outdir/SSR.Rscript";