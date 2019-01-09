#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd;
use Getopt::Long;

my $R="/share/work/biosoft/R/R-v3.0.0/bin/";
my $tgicl="/share/work/biosoft/tgicl/TGICL-2.1/bin/";
my $script_N50="/share/work/fanyc/pipeline/RNAseq_no_ref_pipeline/sub/N50.pl";



my(@sampleID,$outdir,@list,$unigene);
GetOptions(
    'list:s{1,}'=>\@list,
    "o:s"=>\$outdir,
    "fa:s"=>\$unigene,
    "id:s{1,}"=>\@sampleID,
    );
sub usage{
    print qq{
This script plot venn of unigene
usage:
perl $0 -list A01.unigene_id.txt A02.unigene_ID.txt -o /path/to/directory/ -fa unigene.fasta -id A01 A02
-list     the unigene id list,several samples split by space:force
-o      the output directory:froce
-fa     the unigene sequence:froce
-id   the sample ID,several sample split by space:force
};
    exit;
}

if (!@list || !@sampleID || !$unigene || !$outdir)
{
    &usage();
}


my @id=@list;
my @sample=@sampleID;
my @venn;
system `$tgicl/cdbfasta $unigene -o $outdir/unigene.cidx`;
open(OUT, ">$outdir/unigene_statistic.txt");
print OUT "Samples\tNumber of unigenes\tTotal length of unigenes(nt)\tMean length of unigenes(nt)\tN50(nt)\n";

for (my $k=0;$k<=($#id+1);$k++)
{
    my $num;
    my $total_length;
    my $mean_length;
    my $N50;
    my @string;
    my %hash;
    my $seqname;
    if ($k<=$#id)
    {
        `$tgicl/cdbyank $outdir/unigene.cidx < $id[$k] >$outdir/$sample[$k].fasta`;
        $N50=`perl $script_N50 -f $outdir/$sample[$k].fasta -l 0`;
        open(IN,"$outdir/$sample[$k].fasta");
    }
    if ($k==($#id+1))
    {
        $N50=`perl $script_N50 -f $unigene -l 0`;
        open(IN,$unigene);
    }
    chomp($N50);
    @string=split(/ /,$N50);
    while (<IN>)
    {
        chomp;
        if ($_=~/\>/)
        { 
            ++$num;
            $seqname=substr($_,1);
        }
        else
        {
            $hash{$seqname}.=$_;
        }   
    }
    foreach my $key (keys %hash)
    {
        $total_length+=length($hash{$key});
    }
        $mean_length=$total_length/$num;
    if ($k<=$#id)
    {
        print OUT "$sample[$k]\t$num\t$total_length\t$mean_length\t$string[1]\n";
    }
    else
    {
        print OUT "ALL\t$num\t$total_length\t$mean_length\t$string[1]\n";
    } 
    close IN;
}

if ($#id==1)
{
    system "echo '
#!$R/Rscript
library(VennDiagram)
A<-read.table(\"$id[0]\",header=F)
B<-read.table(\"$id[1]\",header=F)
a<-A[,1]
b<-B[,1]
venn.plot <- venn.diagram(
	x = list(
		$sample[0]= a,
        $sample[1]= b
		),
	filename = \"$outdir/Venn.png\",
	lwd = 4,
	fill = c(\"light blue\", \"pink\"),
	alpha = 0.75,
	label.col = \"white\",
	cex = 2,
	fontfamily = \"serif\",
	fontface = \"bold\",
	cat.col = c(\"light blue\", \"pink\"),
	cat.cex = 3,
	cat.fontfamily = \"serif\",
	cat.fontface = \"bold\",
	cat.dist = c(0.03, 0.03),
	cat.pos = c(-20, 14)
	);
    '>$outdir/venn.R";
    system "$R/Rscript $outdir/venn.R";
}
if($#id==2)
{
        system "echo '
#!$R/Rscript
library(VennDiagram) 
A<-read.table(\"$id[0]\",header=F)
B<-read.table(\"$id[1]\",header=F)
C<-read.table(\"$id[2]\",header=F)
a<-A[,1]
b<-B[,1]
c<-C[,1]
venn.plot <- venn.diagram(
	x = list(
	$sample[0]= a,
    $sample[1]= b,
    $sample[2]= c
		),
	filename = \"$outdir/triple_Venn.tiff\",
	col = \"transparent\",
	fill = c(\"red\", \"blue\", \"green\"),
	alpha = 0.5,
	label.col = c(\"darkred\", \"white\", \"darkblue\", \"white\",
	 \"white\", \"white\", \"darkgreen\"),
	cex = 2,
	fontfamily = \"serif\",
	fontface = \"bold\",
	cat.default.pos = \"text\",
	cat.col = c(\"darkred\", \"darkblue\", \"darkgreen\"),
	cat.cex = 2,
	cat.fontfamily = \"serif\",
	cat.dist = c(0.06, 0.06, 0.03),
	cat.pos = 0
	);
    '>$outdir/venn.R";
    system "$R/Rscript $outdir/venn.R";
}
if($#id==3)
{
        system "echo '
#!$R/Rscript
library(VennDiagram)
A<-read.table(\"$id[0]\",header=F)
B<-read.table(\"$id[1]\",header=F)
C<-read.table(\"$id[2]\",header=F)
D<-read.table(\"$id[3]\",header=F)
a<-A[,1]
b<-B[,1]
c<-C[,1]
d<-D[,1]
venn.plot <- venn.diagram(
  x = list(
    $sample[0]= a,
    $sample[1]= b,
    $sample[2]= c,
    $sample[3]= d
  ),
  filename = \"$outdir/Venn.png\",
  col = \"transparent\",
  fill = c(\"cornflowerblue\", \"green\", \"yellow\", \"darkorchid1\"),
  alpha = 0.50,
  label.col = c(\"orange\", \"white\", \"darkorchid4\", \"white\", 
                \"white\", \"white\", \"white\",\"white\",\"darkblue\", \"white\", 
                \"white\", \"white\", \"white\", \"darkgreen\", \"white\"),
  cex = 1.5,
  fontfamily = \"serif\",
  fontface = \"bold\",
  cat.col = c(\"darkblue\", \"darkgreen\", \"orange\", \"darkorchid4\"),
  cat.cex = 1.5,
  cat.pos = 0,
  cat.dist = 0.07,
  cat.fontfamily = \"serif\",
  rotation.degree = 270,
  margin = 0.2
);    
    
    '>$outdir/venn.R";
    system "$R/Rscript $outdir/venn.R";
}
if($#id==4)
{
        system "echo '
#!$R/Rscript
library(VennDiagram)
A<-read.table(\"$id[0]\",header=F)
B<-read.table(\"$id[1]\",header=F)
C<-read.table(\"$id[2]\",header=F)
D<-read.table(\"$id[3]\",header=F)
E<-read.table(\"$id[4]\",header=F)
a<-A[,1]
b<-B[,1]
c<-C[,1]
d<-D[,1]
e<-E[,1]
venn.plot <- venn.diagram(
  x = list(
    $sample[0]= a,
    $sample[1]= b,
    $sample[2]= c,
    $sample[3]= d,
    $sample[4]= e
  ),
  filename = \"$outdir/Venn.png\",
  col = \"black\",
  fill = c(\"dodgerblue\", \"goldenrod1\", \"darkorange1\", \"seagreen3\", \"orchid3\"),
  alpha = 0.50,
  cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
          1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
  cat.col = c(\"dodgerblue\", \"goldenrod1\", \"darkorange1\", \"seagreen3\", \"orchid3\"),
  cat.cex = 1.5,
  cat.fontface = \"bold\",
  margin = 0.05
); 
    '>$outdir/venn.R";
    system "$R/Rscript $outdir/venn.R";
}
