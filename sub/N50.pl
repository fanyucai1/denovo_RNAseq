#! /usr/local/bin/perl5.8.3 -w

###
### Script Name : calcN50.pl
### Description : This script calculates the N50 of a collection of sequences in a fasta file
### Input       : A fasta sequence file and the minimum length of sequence to consider it for the N50 calculation
### Output      : None.  The N50 value as well as the number of sequences >= this value is sent to STDOUT
###

use strict;
use Getopt::Long;
use Bio::SeqIO;

my $fasta_file;
my $filter_length;
GetOptions("fasta_file=s" => \$fasta_file, "length_filter=i" => \$filter_length);
die "Usage: $0 [options]
\t\t-f /path/to/fasta/file [Required]
\t\t-l filter seq less than length (0 if no filter is to be done)\n\n" if !$fasta_file || !defined($filter_length) || $filter_length !~ /\d+/;

###
### Get the length of each sequence:
###
my $in = Bio::SeqIO->new(-file => $fasta_file, -format => 'fasta');
my @seq_lengths = ();
while (my $seq_obj = $in->next_seq){
	next if length($seq_obj->seq) < $filter_length;
	push @seq_lengths, length($seq_obj->seq);
}



###
### Determine the N50:
###
my @vals;
my $sum = 0;
for my $val (@seq_lengths){
    chomp $val;
    $sum += $val;
    push(@vals, $val);
}

my @sorted = sort { $a <=> $b } @vals;

my $runningSum = 0;
my $N50 = 0;
foreach my $v (@sorted)
{
    $runningSum += $v;
    if($runningSum >= ($sum / 2))
    {
	$N50 = $v;
	last;
    }
}


###
### Determine the number of sequences >= than the N50 in length:
###
my $count = 0;
foreach my $v (@sorted)
{
    if($v >= $N50)
    {
	$count++;
    }
}


printf("N50: $N50 Count >= N50: $count\n");


