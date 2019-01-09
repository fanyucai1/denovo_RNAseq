#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin);
use Config::IniFiles;

my $input;
GetOptions(
    "i:s"=>\$input,
           );
sub usage{
    print qq{
This script could analysis 
usage:
perl $0 -i config.ini
options:
-i          input config file,force
Email:fanyucai1\@126.com
2016.3.30
version 1.0
2016.8.29
vesion 2.0
};
    exit;
}
if (!$input )
{
    &usage();
}
sub qsub()
{
	my ($shfile, $queue, $ass_maxproc) = @_ ;
    $queue||="all.q";
    $ass_maxproc||=10;
    my $cmd = "perl /usr/local/bin/qsub-sge.pl --maxproc $ass_maxproc --queue $queue --resource vf=15G --reqsub $shfile --independent" ;
	if (`hostname` !~ /cluster/)
    {
		print "Only cluster could run this process ";
        exit;
	}
    my $flag=system($cmd);
    if($flag !=0)
    {
		die "qsub [$shfile] die with error : $cmd \n";
        exit;
	}
}

my %ini;
tie %ini, 'Config::IniFiles', (-file =>$input);

#mkdir the output directory
foreach my $key (keys %{$ini{dir}})
{
    $ini{dir}{$key}=~ s/\$outdir/$ini{dir}{outdir}/;
    system "mkdir -p $ini{dir}{$key}";
}
#fastqc_Statistics
my ($a,$b,$name);
foreach my $key (keys %{$ini{sample}})
{
    my @array=split(/\|/,$ini{sample}{$key});
    $a.=" $array[0]";
    $b.=" $array[1]";
    $name.=" $key";  
}
system "echo 'perl $Bin/sub/fastqc_Statistics.pl -a $a -b $b -id $name -o $ini{dir}{html_assess}/'>$ini{dir}{shell}/statistics_1.sh";
if (`hostname` !~ /cluster/)
{
    `nohup sh $ini{dir}{shell}/statistics_1.sh&`;
}
else
{
    `nohup perl $ini{sub}{sub} $ini{dir}{shell}/statistics_1.sh``;
}
##########################################trinity assembly
my ($bam);
#combine assembly by trinity
if($ini{config}{assembly} =~"combine")
{
    my ($left_total,$right_total);
    open(TR,">$ini{dir}{shell}/trinity.sh2");
    open(OUT, ">$ini{dir}{shell}/RSEM3.sh");
    my ($type,$last,$venn_list_ID,$sample_ID);
    foreach my $key(keys %{$ini{sample}})
    {
        my @array=split(/\|/,$ini{sample}{$key});
        $left_total.="$array[0],";
        $right_total.="$array[1],";
        print OUT "export PATH=$ini{software}{RSEM}:\$PATH && perl $ini{software}{trinity}/util/align_and_estimate_abundance.pl --transcripts $ini{dir}{outdir}/trinity_out_dir/Trinity.fasta --seqType fq --left $array[0] --right $array[1]  --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference --output_prefix $key --output_dir $ini{dir}{outdir}/rsem_outdir";
        $type.=" $key\.isoforms.results";
        $sample_ID.="$key ";
        $bam.="$ini{dir}{outdir}/rsem_outdir/$key\.bowtie.bam ";
        $last="$ini{dir}{outdir}/rsem_outdir/$key\.isoforms.results";
        $venn_list_ID.="$ini{dir}{unigene}/venn/$key.unigene_ID.txt ";
    }
    chop($left_total);
    chop($right_total);
    print TR "$ini{software}{trinity}/Trinity --no_version_check --inchworm_cpu 20 --bflyHeapSpaceInit 20G --bflyHeapSpaceMax 100G --left $left_total --right $right_total  --CPU $ini{config}{CPU} --output $ini{dir}{outdir}/trinity_out_dir --seqType fq --min_kmer_cov $ini{config}{kmercov} --min_contig_length $ini{config}{contiglenth} --max_memory $ini{config}{memory} && ";
    print TR "perl $ini{software}{trinity}/util/TrinityStats.pl $ini{dir}{outdir}/trinity_out_dir/Trinity.fasta >$ini{dir}{outdir}/trinity_out_dir/contig_Statistics.txt";
    &qsub("$ini{dir}{shell}/trinity2.sh","great.q");
    &qsub("$ini{dir}{shell}/RSEM3.sh");
    
    open(IN, ">$ini{dir}{shell}/unigene4.sh");
    print IN "export PATH=$ini{software}{R}:\$PATH && cd $ini{dir}{outdir}/rsem_outdir/ && ";
    print IN "$ini{software}{trinity}/util/abundance_estimates_to_matrix.pl --est_method RSEM --out_prefix $ini{dir}{DGE}/trans_counts $type && ";
    print IN "perl $Bin/sub/isofrom_get_unigene.pl -i $type -o $ini{dir}{unigene}/ -TPM $ini{config}{TPM} -ID $sample_ID -matrix $ini{dir}{DGE}/trans_counts.counts.matrix && ";
    print IN "mv $ini{dir}{unigene}/matrix.txt $ini{dir}{DGE}/ && ";
    print IN "$ini{software}{trinity}/util/retrieve_sequences_from_fasta.pl $ini{dir}{unigene}/unigene_ID.txt $ini{dir}{outdir}/trinity_out_dir/Trinity.fasta >$ini{dir}{unigene}/unigene.fasta && ";
    print IN "cat $last | cut -f1,3,4 > $ini{dir}{DGE}/length.txt && ";
    print IN "perl $Bin/sub/unigene_venn.pl -list $venn_list_ID -o $ini{dir}{unigene}/venn -fa $ini{dir}{unigene}/unigene.fasta -id $sample_ID";
    close IN;
    &qsub("$ini{dir}{shell}/unigene4.sh");
    open(BAM,">$ini{dir}{shell}/bam_assess.sh");
    print BAM "perl $Bin/sub/bam_stat.pl -bam $bam -o $ini{dir}{assess} -id $sample_ID -bed $ini{dir}{code}/unigene.fasta.transdecoder.bed";  
}
else
{
    my ($left_total,$right_total,$str);
    foreach my $key(keys %{$ini{assembly}})
    {
        open(TR,">$ini{dir}{shell}/trinity2.sh");
        system "mkdir -p $ini{dir}{outdir}/$key/";
        my @name=split(/,/,$ini{assembly}{$key});
        for(my $k=0;$k<=$#name;$k++)
        {
            my @array=split(/\|/,$ini{sample}{$name[$k]});
            $left_total.="$array[0],";
            $right_total.="$array[1],";
        }
        chop($left_total);
        chop($right_total);
        print TR "$ini{software}{trinity}/Trinity --no_version_check --inchworm_cpu 20 --bflyHeapSpaceInit 20G --bflyHeapSpaceMax 100G --left $left_total --right $right_total  --CPU $ini{config}{CPU} --output $ini{dir}{outdir}/$key/trinity_out_dir --seqType fq --min_kmer_cov $ini{config}{kmercov} --min_contig_length $ini{config}{contiglenth} --max_memory $ini{config}{memory}";
        &qsub("$ini{dir}{shell}/trinity.sh2","great.q");
        close;
        $str.="$ini{dir}{outdir}/$key/trinity_out_dir/Trinity.fasta ";
    }
    #use cdhit to delete the same sequence
    system "mkdir -p $ini{dir}{outdir}/cdhit";
    system "cat $str >$ini{dir}{outdir}/cdhit/Trinity_combine.fasta";
    open(CD,">$ini{dir}{shell}/cdhit3.sh");
    print CD "perl $ini{software}{trinity}/util/TrinityStats.pl $ini{dir}{outdir}/cdhit/Trinity_combine.fasta >$ini{dir}{outdir}/trinity_out_dir/contig_Statistics.txt && ";
    print CD "$ini{software}{cdhit} -i $ini{dir}{outdir}/cdhit/Trinity_combine.fasta -o $ini{dir}{outdir}/cdhit/Trinity_cdhit.fasta -c 1 -M 10000 && ";
    print CD "perl $Bin/sub/SolexaQA2qiime.pl -i $ini{dir}{outdir}/cdhit/Trinity_cdhit.fasta -f fa -factor cdhit -qiime $ini{dir}{outdir}/cdhit/cdhit.fasta";
    &qsub("$ini{dir}{shell}/cdhit3.sh");
    #get unigene
    open(UN,">$ini{dir}{shell}/tgicl4.sh");
    system "mkdir -p $ini{dir}{outdir}/tgicl";
    print UN "cd $ini{dir}{outdir}/tgicl && ";
    print UN "export PERL5LIB=$ini{software}{tgicl_lib}:\$PERL5LIB && ";
    print UN "cp $ini{software}{tgicl_config} $ini{dir}{outdir}/tgicl/  && ";
    print UN "$ini{software}{perl} $ini{software}{tgicl}/tgicl -F $ini{dir}{outdir}/cdhit/cdhit.fasta -c $ini{config}{c} -p $ini{config}{p} -l $ini{config}{l} && ";
    print UN "$ini{software}{tgicl}/cdbfasta $ini{dir}{outdir}/cdhit/cdhit.fasta -o $ini{dir}{outdir}/tgicl/cdhit.fasta.cidx && ";
    print UN "$ini{software}{tgicl}/cdbyank $ini{dir}{outdir}/tgicl/cdhit.fasta.cidx < $ini{dir}{outdir}/tgicl/cdhit.fasta.singletons >$ini{dir}{outdir}/tgicl/singletons.fasta && ";
    print UN "cat $ini{dir}{outdir}/tgicl/asm_*/contigs $ini{dir}{outdir}/tgicl/singletons.fasta >$ini{dir}{outdir}/tgicl/tgicl_out.fasta";
    &qsub("$ini{dir}{shell}/tgicl.sh4");
    #prepare the reference for alignment and abundance estimation use kallisto
    open(INDEX,">$ini{dir}{shell}/quant_index5.sh");
    print INDEX "$ini{software}{kallisto}/kallisto index -i $ini{dir}{outdir}/quant/tgicl.idx $ini{dir}{outdir}/tgicl/tgicl_out.fasta";
    &qsub("$ini{dir}{shell}/quant_index5.sh");

    open(MAP,">$ini{dir}{shell}/quant_map6.sh");
    open(OUT, ">$ini{dir}{shell}/unigene7.sh");
    my ($kalliso,@line,$sample_ID,$numb);
    foreach my $key (keys %{$ini{sample}})
    {
        $line[$numb]="$key";
        ++$numb;
        my @array=split(/\|/,$ini{sample}{$key});
        system "mkdir -p $ini{dir}{outdir}/quant/$key/";
        print MAP "$ini{software}{kallisto}/kallisto quant -i $ini{dir}{outdir}/quant/tgicl.idx -o $ini{dir}{outdir}/quant/$key  --pseudobam $array[0] $array[1]| $ini{software}{samtools}/samtools view -Sb - > $ini{dir}{outdir}/quant/$key/$key.bam\n";
        print OUT "cat $ini{dir}{outdir}/quant/$key/abundance.tsv|cut -f1,2,3|sed -e '1c transcript_id\tlength\teffective_length'>$ini{dir}{DGE}/length.txt && ";
        $kalliso.="$ini{dir}{outdir}/quant/$key/abundance.tsv ";
        $sample_ID.="$key ";
        $bam.="$ini{dir}{outdir}/quant/$key/$key.bam ";
    }
    &qsub("$ini{dir}{shell}/quant_map6.sh");
    
    #Counting Numbers of Expressed Transcripts or Genes
    print OUT "export PATH=$ini{software}{R}:\$PATH && perl $ini{software}{trinity}/util/abundance_estimates_to_matrix.pl --est_method kallisto --out_prefix $ini{dir}{DGE}/trans_counts --name_sample_by_basedir $kalliso && ";
    #filter count use TPM and get unigene ID per sample
    my ($linux,$venn_ID,$list_ID);
    for (my $k=0;$k<=$#line;$k++)
    {
        my $i=$k+2;
        print OUT "awk -F\"\\t\" '{if (\$$i>1) print \$1}' $ini{dir}{DGE}/trans_counts.TMM.EXPR.matrix >$ini{dir}{unigene}/venn/$line[$k]\__uigene_ID.txt && ";
        $venn_ID="$ini{dir}{unigene}/venn/$line[$k]\__uigene_ID.txt ";
        $list_ID="$line[$k] ";
        if ($k==0)
        {
            $linux.="\$$i>1";
        }
        else
        {
            $linux.=" ||\$$i>1";
        }
    }
    print OUT "awk -F\"\\t\" '{if ($linux) print \$0}' $ini{dir}{DGE}/trans_counts.TMM.EXPR.matrix >$ini{dir}{DGE}/matrix.txt && ";
    print OUT "cat $ini{dir}{unigene}/venn/*__uigene_ID.txt |sort -u >$ini{dir}{unigene}/unigene_ID.txt && ";
    print OUT "$ini{software}{tgicl}/cdbfasta $ini{dir}{outdir}/tgicl/tgicl_out.fasta -o $ini{dir}{unigene}/tgicl_out.fasta.cidx && ";
    print OUT "$ini{software}{tgicl}/cdbyank $ini{dir}{unigene}/tgicl_out.fasta.cidx < $ini{dir}{unigene}/unigene_ID.txt >$ini{dir}{unigene}/unigene.fasta && ";
    print OUT "perl $Bin/sub/statistics_unigene.pl -fasta $ini{dir}{unigene}/unigene.fasta -o $ini{dir}{unigene}/ && ";
    print OUT "perl $Bin/sub/unigene_venn.pl -list $venn_ID -o $ini{dir}{unigene}/venn/ -f $ini{dir}{unigene}/unigene.fasta -id $list_ID";
    close OUT;
    &qsub("$ini{dir}{shell}/unigene7.sh");
    ########################split_assembly_bam_data assess
    open(BAM,">$ini{dir}{shell}/bam_assess.sh");
    print BAM "perl $Bin/sub/bam_statistics.pl -bam $bam -o $ini{dir}{assess} -id $sample_ID -bed $ini{dir}{code}/unigene.fasta.transdecoder.bed";  
}
#####################DGE analysis
if ($ini{config}{sample_num} >1)
{
    my $num=0;
    foreach my $key ( keys %{$ini{condition}})
    {
        ++$num;
        if ($num==1)
        {
            system "echo '$ini{condition}{$key}\t$key'>$ini{dir}{DGE}/samples.txt";
        }
        else
        {
            system "echo '$ini{condition}{$key}\t$key'>>$ini{dir}{DGE}/samples.txt";
        }
    }
    $num=0;
    my @plot;
    foreach my $key ( keys %{$ini{contrasts}})
    {
        ++$num;
        my @array=split(/,/,$ini{contrasts}{$key});
        if ($num==1)
        {
            system "echo '$array[0]\t$array[1]'>$ini{dir}{DGE}/contrasts.txt";
            system "echo 'perl $Bin/sub/plot_Volcano_MA.pl -i matrix.txt.$array[0]_vs_$array[1].edgeR.DE_results -vs $array[0]_vs_$array[1] -o $ini{dir}{DGE}/edgeR'>$ini{dir}{shell}/Volcano_MA.sh";
        }
        else
        {
            system "echo  '$array[0]\t$array[1]'>>$ini{dir}{DGE}/contrasts.txt";
            system "echo 'perl $Bin/sub/plot_Volcano_MA.pl -i matrix.txt.$array[0]_vs_$array[1].edgeR.DE_results -vs $array[0]_vs_$array[1] -o $ini{dir}{DGE}/edgeR'>>$ini{dir}{shell}/Volcano_MA.sh";
        }
    }
    open(DGE,">$ini{dir}{shell}/DGE.sh");
    print DGE "cd $ini{dir}{DGE} && export PATH=$ini{software}{R}:\$PATH && ";
    print DGE "perl $ini{software}{edgeR_trinity}/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix $ini{dir}{DGE}/matrix.txt --method edgeR --samples_file $ini{dir}{DGE}/samples.txt  --output $ini{dir}{DGE}/edgeR --contrasts $ini{dir}{DGE}/contrasts.txt --dispersion $ini{config}{dispersion} && ";
    print DGE "$ini{software}{edgeR_trinity}/Analysis/DifferentialExpression/run_TMM_normalization_write_FPKM_matrix.pl --lengths $ini{dir}{DGE}/length.txt --matrix $ini{dir}{DGE}/matrix.txt && ";
    print DGE "cd $ini{dir}{DGE}/edgeR && ";
    print DGE "perl $ini{software}{edgeR_trinity}/Analysis/DifferentialExpression/analyze_diff_expr.pl --max_genes_clust 1000 --matrix $ini{dir}{DGE}/matrix.txt.TMM_normalized.FPKM --samples $ini{dir}{DGE}/samples.txt -P $ini{config}{FDR} -C $ini{config}{foldchange}  && "; 
    print DGE "export PATH=$ini{software}{R3_2}:\$PATH && ";
    print DGE "sed -i s/ColSideColors=sampleAnnotations\\)/ColSideColors=sampleAnnotations,labRow=F\\)/g *.matrix.R && ";
    print DGE "sed -i s/pdf/png/g  *.matrix.R && ";
    print DGE "$ini{software}{R}/Rscript  *.matrix.R && ";
    print DGE "$ini{software}{perl} $Bin/sub/FPKM_boxplot_density_plot.pl -i $ini{dir}{DGE}/matrix.txt.TMM_normalized.FPKM -o $ini{dir}{DGE}/ && ";
    print DGE "sh $ini{dir}{shell}/Volcano_MA.sh";
    &qsub("$ini{dir}{shell}/DGE.sh");
}

#snp analysis
open(SNPI,">$ini{dir}{shell}/SNP_index.sh");
print SNPI "cd $ini{dir}{unigene}/ && ";
print SNPI "$ini{software}{java}  -XX:ParallelGCThreads=5 -jar $ini{software}{picard} CreateSequenceDictionary REFERENCE=unigene.fasta OUTPUT=unigene.dict && ";
print SNPI "$ini{software}{STAR}/STAR  --runMode genomeGenerate --runThreadN 10 --genomeDir $ini{dir}{SNP}/ --genomeFastaFiles $ini{dir}{unigene}/unigene.fasta --limitGenomeGenerateRAM 70000000000 && ";
print SNPI "$ini{software}{samtools}/samtools faidx unigene.fasta";
&qsub("$ini{dir}{shell}/SNP_index.sh","great.q");

open(SNP, ">$ini{dir}{shell}/SNP_map.sh");
my $linux=$ini{config}{readlength}-1;
my $string;
foreach my $key(keys %{$ini{sample}})
{
    my @array=split(/\|/,$ini{sample}{$key});
    system "mkdir -p $ini{dir}{SNP}/$key";
    print SNP "cd $ini{dir}{SNP}/$key/ && ";
    print SNP "$ini{software}{STAR}/STAR --genomeDir $ini{dir}{SNP} --runThreadN 5 --sjdbOverhang $linux --readFilesIn $array[0] $array[1] --outSAMtype BAM SortedByCoordinate --alignSoftClipAtReferenceEnds No --limitBAMsortRAM  70000000000 --twopassMode Basic && "; 
    #Add read groups, sort, mark duplicates, and create index
    print SNP "$ini{software}{java} -Xmx10g  -XX:ParallelGCThreads=5 -jar $ini{software}{picard} AddOrReplaceReadGroups INPUT=$ini{dir}{SNP}/$key/Aligned.sortedByCoord.out.bam OUTPUT=$ini{dir}{SNP}/$key/addRG.bam RGID=$key RGLB=$key RGPL=illumina RGPU=$key RGSM=$key SORT_ORDER=coordinate CREATE_INDEX=true && ";
    print SNP "$ini{software}{java} -Xmx10g  -XX:ParallelGCThreads=5 -jar $ini{software}{picard} MarkDuplicates CREATE_INDEX=true INPUT=$ini{dir}{SNP}/$key/addRG.bam OUTPUT=$ini{dir}{SNP}/$key/markdup.bam METRICS_FILE=$ini{dir}{SNP}/$key/metrics.txt && ";
    #Split'N'Trim and reassign mapping qualities
    print SNP "$ini{software}{java} -Xmx10g  -XX:ParallelGCThreads=5 -jar $ini{software}{gatk} -T SplitNCigarReads -R $ini{dir}{unigene}/unigene.fasta  -I $ini{dir}{SNP}/$key/markdup.bam -o $ini{dir}{SNP}/$key/split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS";
    $string.="$ini{dir}{SNP}/$key/split.bam ";
}
&qsub("$ini{dir}{shell}/SNP_map.sh");
close SNP;
open(CALL, ">$ini{dir}{shell}/call_snp.sh");
    #Variant calling
print CALL "perl $Bin/sub/call_snp.pl -chr $ini{dir}{SNP}/chrName.txt -bam $string -ref $ini{dir}{unigene}/unigene.fasta -num 2000 -o $ini{dir}{SNP}/ && ";
    #Variant filtering
print CALL "$ini{software}{java} -Xmx40g  -Djava.io.tmpdir=/share/work/tmp -XX:ParallelGCThreads=5 -jar $ini{software}{gatk} -T VariantFiltration -R $ini{dir}{unigene}/unigene.fasta -V $ini{dir}{SNP}/raw.vcf -window 35 -cluster 3 -filterName FS -filter \"FS > 30.0\" -filterName QD -filter \"QD < 2.0\" -o $ini{dir}{SNP}/snp_indels.vcf && ";
    #Extract the SNPs from the call set
print CALL "$ini{software}{java} -Xmx40g  -Djava.io.tmpdir=/share/work/tmp -XX:ParallelGCThreads=5 -jar $ini{software}{gatk} -T SelectVariants -R $ini{dir}{unigene}/unigene.fasta -V $ini{dir}{SNP}/snp_indels.vcf -selectType SNP -o $ini{dir}{SNP}/snp.vcf && ";
    #Extract the Indels from the call set
print CALL "$ini{software}{java} -Xmx40g  -Djava.io.tmpdir=/share/work/tmp -XX:ParallelGCThreads=5 -jar $ini{software}{gatk} -T SelectVariants -R $ini{dir}{unigene}/unigene.fasta -V $ini{dir}{SNP}/snp_indels.vcf -selectType INDEL -o $ini{dir}{SNP}/indels.vcf && ";
    #Extract genotypes in variant sites shared by all individuals
print CALL "cd $ini{dir}{SNP}/ && ";
print CALL "$Bin/sub/getgenosfromvcf.py $ini{dir}{SNP}/snp_indels.vcf $ini{dir}{SNP}/Genotypes.txt rows 20 && ";
    #create a new file containing only the variant sites for which we have genotype information for all individuals
print CALL "cat $ini{dir}{SNP}/Genotypes.txt | grep -v '\.' > $ini{dir}{SNP}/genotypes_shared_by_all.txt";
close CALL;
&qsub("$ini{dir}{shell}/call_snp.sh");

#CDS analysis
system "echo 'cd $ini{dir}{code} && $ini{software}{TransDecoder}/TransDecoder.LongOrfs -t $ini{dir}{unigene}/unigene.fasta'>$ini{dir}{shell}/code.sh";
system "echo '$ini{software}{hmmscan} --cpu 8 --domtblout $ini{dir}{code}/pfam.domtblout $ini{annotation}{Pfam} $ini{dir}{code}/unigene.fasta.transdecoder_dir/longest_orfs.pep'>>$ini{dir}{shell}/code.sh";
system "echo '$ini{software}{blastp} -query $ini{dir}{code}/unigene.fasta.transdecoder_dir/longest_orfs.pep  -db $ini{annotation}{Swissprot} -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 10 > $ini{dir}{code}/blastp.outfmt6'>>$ini{dir}{shell}/code.sh";
system "echo '$ini{software}{TransDecoder}/TransDecoder.Predict  -t $ini{dir}{unigene}/unigene.fasta --retain_pfam_hits $ini{dir}{code}/pfam.domtblout --retain_blastp_hits $ini{dir}{code}/blastp.outfmt6'>>$ini{dir}{shell}/code.sh";
$linux=`wc -l $ini{dir}{shell}/code.sh`;
chomp($linux);
system "perl $ini{qsub}{qsub} -line $linux";
&qsub("$ini{dir}{shell}/bam_assess.sh");

#SSR analysis
system "echo 'perl $ini{software}{MISA} $ini{dir}{unigene}/unigene.fasta'>$ini{dir}{shell}/SSR.sh";
system "echo 'mv $ini{dir}{unigene}/*statistics $ini{dir}{unigene}/*misa $ini{dir}{SSR}/'>>$ini{dir}{shell}/SSR.sh";
system "echo 'cp  $ini{dir}{unigene}/unigene.fasta $ini{dir}{SSR}/'>>$ini{dir}{shell}/SSR.sh";
system "echo 'cd $ini{dir}{SSR}/'>>$ini{dir}{shell}/SSR.sh";
system "echo '$ini{software}{perl} $Bin/sub/p3_in.pl unigene.fasta.misa unigene.fasta unigene.fasta.p3in'>>$ini{dir}{shell}/SSR.sh";
system "echo '$ini{software}{primer3}  -default_version=1 -io_version=3 <unigene.fasta.p3in > unigene.fasta.p3out '>>$ini{dir}{shell}/SSR.sh";
system "echo '$ini{software}{perl} $Bin/sub/p3_out.pl unigene.fasta.p3out unigene.fasta.misa'>>$ini{dir}{shell}/SSR.sh";
system "echo '$ini{software}{perl} $Bin/sub/plot_SSR.pl -i $ini{dir}{SSR}/unigene.fasta.statistics -o $ini{dir}{SSR}/'>>$ini{dir}{shell}/SSR.sh";
$linux=`wc -l $ini{dir}{shell}/SSR.sh`;
chomp($linux);
system "perl $ini{qsub}{qsub} -line $linux";

#statistics_unigene
system "echo 'nohup perl $Bin/sub/statistics_unigene.pl -fasta $ini{dir}{unigene}/unigene.fasta -o $ini{dir}{unigene}/ &'>>$ini{dir}{shell}/bin.sh";

#annotation
open(OUTT, ">$ini{dir}{annotation}/anno_config.txt");
my $data;
print OUTT "mRNA\t$ini{dir}{unigene}/unigene.fasta\n";
foreach my $key(keys %{$ini{annotation}})
{
    if ($key =~/Cog/i)
    {
        $data="--cog";    
    }
    if($key=~/Kog/i)
    {
        $data="--kog";
    }
    print OUTT "$key\t$ini{annotation}{$key}\n";
}
system "$ini{software}{perl} $ini{software}{anno_script} --nr --swissprot $data --kegg --pfam --GO --trembl --cfg $ini{dir}{annotation}/anno_config.txt -od $ini{dir}{annotation}/";
close OUTT;
#enrichment analysis
if ($ini{config}{sample_num}>1)
{
    open(OUT,">$ini{dir}{shell}/enrichment.sh");
    foreach my $key(keys %{$ini{contrasts}})
    {
        my @array=split(/,/,$ini{contrasts}{$key});
        system "mkdir -p $ini{dir}{enrichment}/$array[0]_vs_$array[1]/";
        print OUT "cat $ini{dir}{DGE}/edgeR/matrix.txt.$array[0]_vs_$array[1].edgeR.DE_results.P$ini{config}{FDR}_C$ini{config}{foldchange}.$array[0]-UP.subset| awk  \'{print \$1\"\\t0\\t0\\t0\\t\"\$5\"\\t\"\$2\"\\tup\"}\'|sed \"1c id   0   0   0   FDR log2FC  regulated\">$ini{dir}{DGE}/edgeR/$array[0]_vs_$array[1]_enrichment_input.txt\n";
        print OUT "awk  \'{print \$1\"\\t0\\t0\\t0\\t0\\t\"\$2\"\\tdown\"}\' $ini{dir}{DGE}/edgeR/matrix.txt.$array[0]_vs_$array[1].edgeR.DE_results.P$ini{config}{FDR}_C$ini{config}{foldchange}.$array[1]-UP.subset|sed 1d>>$ini{dir}{DGE}/edgeR/$array[0]_vs_$array[1]_enrichment_input.txt\n";
        print OUT "$ini{software}{perl} $ini{software}{enrichment_1} -d $ini{dir}{DGE}/edgeR/$array[0]_vs_$array[1]_enrichment_input.txt -k $array[0]_vs_$array[1] -i $ini{dir}{annotation}/Result -o $ini{dir}{enrichment}/$array[0]_vs_$array[1]/ -map /share/work/database/blastdb/db_file.cfg\n";
        print OUT "$ini{software}{perl} $ini{software}{enrichment_2} -enrich_file $ini{dir}{enrichment}/$array[0]_vs_$array[1]/pathway/kegg_enrichment/$array[0]_vs_$array[1].KEGG.enrich.stat.xls -od $ini{dir}{enrichment}/$array[0]_vs_$array[1]/pathway/kegg_enrichment/ -key $array[0]_vs_$array[1]\n";
        print OUT "$ini{software}{perl} $ini{software}{enrichment_3} -All_GO $ini{dir}{enrichment}/$array[0]_vs_$array[1]/go_enrichment/$array[0]_vs_$array[1].GO.list.txt -DEG_list $ini{dir}{DGE}/edgeR/$array[0]_vs_$array[1]_enrichment_input.txt -od $ini{dir}{enrichment}/$array[0]_vs_$array[1]/go_enrichment/ -key $array[0]_vs_$array[1]\n";
        print OUT "$ini{software}{perl} $ini{software}{enrichment_4} -indir $ini{dir}{enrichment}/$array[0]_vs_$array[1]/go_enrichment/ -link $ini{software}{goname}\n";
        print OUT "$ini{software}{R}/Rscript $ini{software}{enrichment_5} $ini{dir}{enrichment}/$array[0]_vs_$array[1]/go_enrichment/  $array[0]_vs_$array[1]\n";
        print OUT "$ini{software}{perl} $ini{software}{enrichment_6} --ipf $ini{dir}{enrichment}/$array[0]_vs_$array[1]/pathway/kegg_enrichment/$array[0]_vs_$array[1].KEGG.xls --opd $ini{dir}{enrichment}/$array[0]_vs_$array[1]/pathway/kegg_enrichment/ --prf $array[0]_vs_$array[1].KEGG\n";
    }
    close OUT;
    $linux=`wc -l $ini{dir}{shell}/enrichment.sh`;
    chomp($linux);
    system "perl $ini{qsub}{qsub} -line $linux";
}
################copy the result to html
system "echo 'cp $Bin/sub/no_ref_RNAseq.png $ini{dir}{html_pic}/'>$ini{dir}{shell}/html.sh";
system "echo 'cp -r $Bin/sub/method_without_reference.pdf   $ini{dir}{html_pic}/'>>$ini{dir}{shell}/html.sh";
#sequencing assess
system "echo 'cp -r $ini{dir}{assess}/*.pdf $ini{dir}{html_assess}/'>>$ini{dir}{shell}/html.sh";
system "echo 'cp -r $ini{dir}{assess}/*.png $ini{dir}{html_assess}/'>>$ini{dir}{shell}/html.sh";

#unigene
system "echo 'cp $ini{dir}{unigene}/unigene.fasta $ini{dir}{html_unigene}'>>$ini{dir}{shell}/html.sh";
system "echo 'cp $ini{dir}{unigene}/Unigene_length.pdf $ini{dir}{unigene}/Unigene_length.png $ini{dir}{html_unigene}'>>$ini{dir}{shell}/html.sh";
if ($ini{config}{sample_num} <=5 && $ini{config}{sample_num}>1)
{
    system "echo 'cp $ini{dir}{unigene}/venn/*png $ini{dir}{html_unigene}'>>$ini{dir}{shell}/html.sh";
}
if ($ini{config}{assembly} =~"combine")
{
    system "echo 'cp $ini{dir}{outdir}/trinity_out_dir/contig_Statistics.txt $ini{dir}{html_unigene}'>>$ini{dir}{shell}/html.sh";
}
system "echo 'cp $ini{dir}{unigene}/venn/*fasta $ini{dir}{html_unigene}'>>$ini{dir}{shell}/html.sh";

#SSR
system "echo 'cp $ini{dir}{SSR}/primer.txt $ini{dir}{html_SSR}'>>$ini{dir}{shell}/html.sh";
system "echo 'cp $ini{dir}{SSR}/*png $ini{dir}{html_SSR}'>>$ini{dir}{shell}/html.sh";
system "echo 'cp $ini{dir}{SSR}/*statistics $ini{dir}{html_SSR}'>>$ini{dir}{shell}/html.sh";
system "echo 'cp $ini{dir}{SSR}/*misa $ini{dir}{html_SSR}'>>$ini{dir}{shell}/html.sh";

#code
system "echo 'cp -r $ini{dir}{code}/* $ini{dir}{html_code}/'>>$ini{dir}{shell}/html.sh";

#SNP
system "echo 'cp $ini{dir}{SNP}/raw.vcf $ini{dir}{html_SNP}/'>>$ini{dir}{shell}/html.sh";
system "echo 'cp $ini{dir}{SNP}/snp_indels.vcf $ini{dir}{html_SNP}/'>>$ini{dir}{shell}/html.sh";
system "echo 'cp $ini{dir}{SNP}/snp.vcf $ini{dir}{html_SNP}/'>>$ini{dir}{shell}/html.sh";
system "echo 'cp $ini{dir}{SNP}/indels.vcf $ini{dir}{html_SNP}/'>>$ini{dir}{shell}/html.sh";
system "echo 'cp $ini{dir}{SNP}/Genotypes.txt $ini{dir}{html_SNP}/'>>$ini{dir}{shell}/html.sh";

if ($ini{config}{sample_num}>1)
{   #DGE
    system "echo 'cp $ini{dir}{DGE}/matrix.txt $ini{dir}{html_DGE}'>>$ini{dir}{shell}/html.sh";
    system "echo 'cp $ini{dir}{DGE}/matrix.txt.TMM_normalized.FPKM $ini{dir}{html_DGE}'>>$ini{dir}{shell}/html.sh";
    system "echo 'cp $ini{dir}{DGE}/edgeR/*png $ini{dir}{html_DGE}'>>$ini{dir}{shell}/html.sh";
    system "echo 'cp $ini{dir}{DGE}/edgeR/*xls $ini{dir}{html_DGE}'>>$ini{dir}{shell}/html.sh";
    system "echo 'cp $ini{dir}{DGE}/*png $ini{dir}{html_DGE}'>>$ini{dir}{shell}/html.sh";
    system "echo 'cp $ini{dir}{DGE}/*pdf $ini{dir}{html_DGE}'>>$ini{dir}{shell}/html.sh";
    #enrichment
    system "echo 'cp -r $ini{dir}{enrichment}/* $ini{dir}{html_enrichment}'>>$ini{dir}{shell}/html.sh";
    system "echo 'cp -r $Bin/sub/enrichment_readme.pdf   $ini{dir}{html_enrichment}/'>>$ini{dir}{shell}/html.sh";
}
else
{
    system "echo 'rm $ini{dir}{html_enrichment}'>>$ini{dir}{shell}/html.sh";
    system "echo 'rm $ini{dir}{html_DGE}'>>$ini{dir}{shell}/html.sh";
}

#annotation
system "echo 'cp -r $ini{dir}{annotation}/Result/* $ini{dir}{html_anno}'>>$ini{dir}{shell}/html.sh";

##get the result
system "sh $ini{dir}{shell}/html.sh";
system "cd $ini{dir}{html} && $ini{software}{perl} $Bin/sub/xml_report.pl -i $input -o $ini{dir}{html}'>>$ini{dir}{shell}/bin.sh";
system "cd $ini{dir}{html} && $ini{software}{python} $ini{software}{xml2html} -i $ini{dir}{html}/output.xml -o $ini{dir}{html}/'>>$ini{dir}{shell}/bin.sh";

