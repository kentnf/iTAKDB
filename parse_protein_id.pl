#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;

my $usage = qq'
usage: perl parse_protein_id.pl [input_seq] [num]

* default sperater is "|";
* example: if The gene is is AAA|BBB|gene_id, the number must be 3;

';

my $file = shift || die $usage;
my $num = shift || 1;

if ($num < 1) { die "Error in num: $num\n"; }

my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$file);
while(my $inseq = $in->next_seq)
{
        print $inseq->id,"\n";
        print $inseq->desc,"\n";

	my @ids = split(/\|/, $inseq->id);
	my $id = $ids[$num-1];
	print ">".$id."\n".$inseq->seq."\n";

}
