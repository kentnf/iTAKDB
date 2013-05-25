#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;
use IO::File;

my $usage = qq'
usage: $0  previous_TFs_seq  previous_PKs_seq  new_TF_seq  new_PK_seq

';

my $preTFs = shift || die $usage;
my $prePKs = shift || die $usage;
my $newTFs = shift || die $usage;
my $newPKs = shift || die $usage;

my %seqHash;

my $in1 = Bio::SeqIO->new(-format=>'fasta', -file=>$preTFs);
while(my $inseq = $in1->next_seq) { $seqHash{$inseq->id} = 1; }

my $in2 = Bio::SeqIO->new(-format=>'fasta', -file=>$prePKs);
while(my $inseq = $in2->next_seq) { $seqHash{$inseq->id} = 1; }

my $output_seq = "Prepared_Protein.fa";
my $out = IO::File->new(">".$output_seq) || die "Can not open output $output_seq\n";

my $in3 = Bio::SeqIO->new(-format=>'fasta', -file=>$newTFs);
while(my $inseq = $in3->next_seq)
{
	unless ( defined $seqHash{$inseq->id} )
	{
		print $out ">".$inseq->id." ".$inseq->desc."\n".$inseq->seq."\n";
		$seqHash{$inseq->id} = 1;
	}
}

my $in4 = Bio::SeqIO->new(-format=>'fasta', -file=>$newPKs);
while(my $inseq = $in4->next_seq)
{
	unless ( defined $seqHash{$inseq->id} )
	{
		print $out ">".$inseq->id." ".$inseq->desc."\n".$inseq->seq."\n";
		$seqHash{$inseq->id} = 1;
	}
}

$out->close;
