#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Cwd;

my ($chr_name,@bam,$unigene,$thread,$num,$outdir);
$num ||=2000;
GetOptions(
    "chr:s"=>\$chr_name,
    "bam:s{1,}"=>\@bam,
    "ref:s"=>\$unigene,
    "t:s"=>\$thread,
    "num:s"=>\$num,
    "o:s"=>\$outdir,
           );
my $gatk="/share/work/biosoft/GATK/3.4.0/GenomeAnalysisTK.jar";
my $java="/share/work/biosoft/java/jre1.8.0_73/bin/java";
my $qsub="/usr/local/bin/qsub-sge.pl";

sub usage{
    print qq{
This script will call snp by split chr.
usage:
perl $0 -chr chr_name.txt -bam sample1.bam,sample2.bam -ref unigene.fasta -num 2000 -o /path/to/outdir/
-chr            the chr name
-bam            the bam file
-ref            the reference file fasta
-num            the split number per process
-o              the outputdirectory
    };
    exit;
}

if (!$chr_name || !@bam || !$unigene || !$outdir)
{
    &usage();
}
system "mkdir -p $outdir/split_snp";

my @bam_array=@bam;
my $string;
for (my $j=0;$j<=$#bam_array;$j++)
{
    $string.=" -I $bam_array[$j] ";
}

open(NAME,$chr_name);
my $line=0;
my $k=1;
while (<NAME>)
{
    chomp;
    ++$line;
    if ($line<=$num)
    {
        system "echo '$_'>>$outdir/split_snp/chr_name_list_$k.list";
    }
    else
    {
        ++$k;
        $line=0;
        system "echo '$_'>$outdir/split_snp/chr_name_list_$k.list";
    }
}
my $sh=1;

for (my $i=1;$i<=$k;$i++)
{
    my $threshold=$sh*10;
    if ($i==1)
    {
        system "echo '$java -Xmx10g  -Djava.io.tmpdir=/share/work/tmp -XX:ParallelGCThreads=10 -jar $gatk -T HaplotypeCaller -nct 10 -R $unigene $string -dontUseSoftClippedBases -L $outdir/split_snp/chr_name_list_$i.list -stand_call_conf 20.0 -stand_emit_conf 20.0 -o $outdir/split_snp/raw_vcf_$i.vcf'>$outdir/split_snp/call_snp_split_$sh.sh";
    }
    elsif ($i>1 && $i<=$threshold)
    {
        system "echo '$java -Xmx10g  -Djava.io.tmpdir=/share/work/tmp -XX:ParallelGCThreads=10 -jar $gatk  -T HaplotypeCaller -nct 10 -R $unigene $string -dontUseSoftClippedBases -L $outdir/split_snp/chr_name_list_$i.list -stand_call_conf 20.0 -stand_emit_conf 20.0 -o $outdir/split_snp/raw_vcf_$i.vcf'>>$outdir/split_snp/call_snp_split_$sh.sh";
    }
    else
    {
       ++$sh;
       system "echo '$java -Xmx10g  -Djava.io.tmpdir=/share/work/tmp -XX:ParallelGCThreads=10 -jar $gatk  -T HaplotypeCaller -nct 10 -R $unigene $string -dontUseSoftClippedBases -L $outdir/split_snp/chr_name_list_$i.list -stand_call_conf 20.0 -stand_emit_conf 20.0 -o $outdir/split_snp/raw_vcf_$i.vcf'>$outdir/split_snp/call_snp_split_$sh.sh";
    }
}

for(my $i=1;$i<=$sh;$i++)
{
        
    `sh $outdir/split_snp/call_snp_split_$i.sh`;
}
#combine vcf
my @ahgvcf;
for(my $i=0;$i<$k;$i++)
{
    my $j=$i+1;
    $ahgvcf[$i]="$outdir/split_snp/raw_vcf_$j.vcf";
}
system "echo 'cat $ahgvcf[0] |grep \"#\" >$outdir/split_snp/head && cat  @ahgvcf |grep -v \"#\" >$outdir/split_snp/body && cat $outdir/split_snp/head $outdir/split_snp/body >$outdir/split_snp/raw.vcf && rm $outdir/split_snp/body $outdir/split_snp/head'>$outdir/split_snp/combine.sh";
system "sh $outdir/split_snp/combine.sh";
system "cp $outdir/split_snp/raw.vcf  $outdir/";