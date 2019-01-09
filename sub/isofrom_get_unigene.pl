#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd;
use Getopt::Long;

my(@input,$outdir,$fasta,@ID,$matrix);
my $threshold||=1;
GetOptions(
    "i:s{1,}"=>\@input,
    "o:s"=>\$outdir,
    "fa:s"=>\$fasta,
    "id:s{1,}"=>\@ID,
    "TPM:s"=>\$threshold,
    "matrix:s"=>\$matrix,
           );

sub usage{
    print qq{
This script you input the isoforms rseult from trinity,you will get unigene and matrix.
perl $0 -i A01.isoforms.results A02.isoforms.results A03.isoforms.results -o /path/to/outdir/ -fa Trinity.fasta -id A01 A02 A03

-i          the RSEM rseult from trinity,several files split by space
-o          the output directory
-ID         the sample ID corrspnoding the input,split by space
-TPM        default,remove if transcripts per million (TPM) <1
-matrix     the estimated RNA-Seq fragment counts (raw counts) from trinity/util/abundance_estimates_to_matrix.pl
Email:fanyucai1\@126.com
2016.8.29
    };
    exit;
}
system "mkdir -p $outdir/venn";
my @array=@input;
my @sampleID=@ID;
for(my $i=0;$i<=$#array;$i++)
{
    open(IN,$array[$i]);
    open(OUT, ">$outdir/$sampleID[$i]\.isoforms.results");
    print OUT "transcript_id	gene_id	length	effective_length	expected_count	TPM	FPKM	IsoPct\n";
    my %hash1; my %hash2;
    my $num;
    while (<IN>)
    {
        chomp;
        ++$num;
        if ($num>1)
        {
            my @temp=split(/\t/,$_);
            if ( $temp[5]>$threshold)
            {
                if ( !$hash1{$temp[1]}||$hash1{$temp[1]}<$temp[2] )
                {
                    $hash1{$temp[1]}=$temp[2];
                    $hash2{$temp[1]}=$_;
                }
            }   
        }
    }
    foreach my $key (keys %hash2)
    {
        print OUT $hash2{$key},"\n";
    }
    close IN;
    close OUT;
    system "cat $outdir/$sampleID[$i]\.isoforms.results|awk -F\"\t\" \'{print \$1}\'|sed \'1d\'>$outdir/venn/$sampleID[$i].unigene_ID.txt";
    if ($i==0)
    {
        system "cat $outdir/venn/$sampleID[$i].unigene_ID.txt >$outdir/all_unigene_ID.txt";
    }
    else
    {
        system "cat $outdir/venn/$sampleID[$i].unigene_ID.txt >>$outdir/all_unigene_ID.txt";
    }
}
system "cat $outdir/all_unigene_ID.txt|sort -u >$outdir/unigene_ID.txt";
system "rm $outdir/all_unigene_ID.txt";
open(IN,"$outdir/unigene_ID.txt");
open(INN, $matrix);
my %unigene;
while (<IN>) {
    chomp;
    $unigene{$_}=1;  
}
my $num=0;
while (<INN>) {
    chomp;
    ++$num;
    my @arra=split;
    if ($num==1)
    {
        system "echo '$_'>$outdir/matrix.txt";
    }
    if ($unigene{$arra[0]})
    {
        system "echo '$_'>>$outdir/matrix.txt";
    }  
}





