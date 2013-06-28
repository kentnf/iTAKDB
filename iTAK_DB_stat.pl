#!/usr/bin/perl

=head1

 update the iTAK database using itak results

=cut
use strict;
use warnings;
use FindBin;
use IO::File;
use Bio::SeqIO;
use File::Basename;
use Getopt::Long;

my $usage = qq'
usage : $0 -i list

  * the format of list file: 
  1) protein sequence
  2) cotyledon ( monocotyledon, dicotyledon, non-angiosperms )
  3) species
  4) transcript sequences
  5) transcript gene file

  * the format of transcript gene file
  1) protein id
  2) transcript id
  3) gene id

';

my ($help, $list);

GetOptions(
	"h"	=> \$help,
	"i=s"	=> \$list,
);

die $usage if $help;
die $usage unless $list;

#================================================================
# check the list file
#================================================================
my %species; my ($has_dicotyledon, $has_monocotyledon) = (0, 0);

my $cfh = IO::File->new($list) || die "Can not open list file $list $!\n";
while(<$cfh>)
{
	chomp; 
	if ($_ =~ m/^#/) { next; }
	my ($pep_file, $cotyledon, $species, $cds_file, $trans_gene_file, $pk_aln, $pk_cat, $pk_seq, $tf_aln, $tf_cat, $tf_seq);
	($pep_file, $cotyledon, $species, $cds_file, $trans_gene_file) = split(/\t/, $_);

	# check cotyledon
	if ($cotyledon eq "monocotyledon" ) 		{ $has_monocotyledon = 1; }
	elsif ( $cotyledon eq "dicotyledon" ) 		{ $has_dicotyledon = 1; }
	elsif ( $cotyledon eq "non-angiosperms" ) 	{}
	else { die "Error at cotyledon info for $_\n"; }

	# species to hash and check the repeat species
	if (defined $species{$species}) { die "Error, repeat species $species\n"; } else { $species{$species} = 1;}

	# set iTAK output files
	# output folder is: TAIR9_protein_output
	# output files are: TAIR9_protein_pkaln, TAIR9_protein_pkcat, TAIR9_protein_pkseq
	#                   TAIR9_protein_tf_align, TAIR9_protein_tf_family, TAIR9_protein_tf_seq
	my $folder = $pep_file."_output";
	my ($fname,$fpath,$fsuffix) = fileparse($pep_file);
	$pk_aln = $folder."/".$fname."_pkaln";
	$pk_cat = $folder."/".$fname."_pkcat";
	$pk_seq = $folder."/".$fname."_pkseq";
	$tf_aln = $folder."/".$fname."_tf_align";
	$tf_cat = $folder."/".$fname."_tf_family";
	$tf_seq = $folder."/".$fname."_tf_seq";

	# check files of iTAK results
	my $run_itak = 0;
	my @files = ($pk_aln, $pk_cat, $pk_seq, $tf_aln, $tf_cat, $tf_seq, $pep_file, $cds_file, $trans_gene_file);
	foreach my $file ( @files ) {
		unless (-s $file) { print "Error! File $file do not exist or has content!\n"; $run_itak = 1;}
	}

	# run itak if the iTAK result is not exist
	if ($run_itak)
	{
		my $itak_cmd = "iTAK.pl -i $pep_file";
		print $itak_cmd."\n";
		#system($itak_cmd) && die "Error in iTAK command: $itak_cmd\n";
	}
}
$cfh->close;

if ($has_dicotyledon) 	{ unless (defined $species{"Arabidopsis"}) { die "Error, no Arabidopsis\n"; } }
if ($has_monocotyledon)	{ unless (defined $species{"Rice"}) { die "Error, no Rice\n"; } }

# check the cluster program
my $clustalw2_program = ${FindBin::RealBin}."/bin/clustalw2";
unless ( -s $clustalw2_program ) { die "Error, can not locate clustalw2 program in $clustalw2_program\n"; }

#================================================================
# list file into to hash					
#================================================================
# init hash for clustering : key1: cyto; key2: species; key3: TF or PK; value: genes (tab delimit)
my %itak_obj;

# init hash for sequences
my %protein_seq;	# key: protein_id,    value: sequence
my %transcript_seq;	# key; transcript_id, value: sequence

# init hash for databases
my %protein_domain;	# key: protein_id, value: domain alignment
my %family_sum;		# key: family_Name, value: number of genes

my %protein_family;	# key: protein_id, value: family
my %protein_species;	# key: protein_id, value: species
my %species_category;	# key: species, value: cotyledon category
my %trans_gene;		# key: transcript_id, value: gene_id
my %trans_protein;	# key: transcript_id, value: protein_id
my %protein_trans;	# key: protein_id, value: transcript_id

my $fh = IO::File->new($list) || die "Can not open list file $list $!\n";
while(<$fh>)
{
	chomp;
	if ($_ =~ m/^#/) { next; }

	# set input files
	my ($pep_file, $cotyledon, $species, $cds_file, $trans_gene_file, $pk_aln, $pk_cat, $pk_seq, $tf_aln, $tf_cat, $tf_seq);
	($pep_file, $cotyledon, $species, $cds_file, $trans_gene_file) = split(/\t/, $_);

	# set iTAK output file
	# output folder is: TAIR9_protein_output
	# output files are: TAIR9_protein_pkaln, TAIR9_protein_pkcat, TAIR9_protein_pkseq
	#                   TAIR9_protein_tf_align, TAIR9_protein_tf_family, TAIR9_protein_tf_seq
	my $folder = $pep_file."_output";
        my ($fname,$fpath,$fsuffix) = fileparse($pep_file);
        $pk_aln = $folder."/".$fname."_pkaln";
        $pk_cat = $folder."/".$fname."_pkcat";
        $pk_seq = $folder."/".$fname."_pkseq";
        $tf_aln = $folder."/".$fname."_tf_align";
        $tf_cat = $folder."/".$fname."_tf_family";
        $tf_seq = $folder."/".$fname."_tf_seq";

	# classification to hash
	my $family;
	my $pcat = IO::File->new($pk_cat) || die "Can not open pk cat file $pk_cat $!\n";
	my $tcat = IO::File->new($tf_cat) || die "Can not open tf cat file $tf_cat $!\n";
	while(<$pcat>) 
	{
		chomp;
		my @a = split(/\t/, $_);
		$family = $a[1];

		if (defined $itak_obj{$cotyledon}{$species}{$family}) {
			$itak_obj{$cotyledon}{$species}{$family}.="\t".$a[0];
		} else {
			$itak_obj{$cotyledon}{$species}{$family} = $a[0];
		}

		if (defined $protein_family{$a[0]}) {
                        $protein_family{$a[0]}.="\t".$a[1];
		} else {
			$protein_family{$a[0]} = $a[1];
		}

		$protein_species{$a[0]} = $species;
	}
	while(<$tcat>)
	{
		chomp;
		my @a = split(/\t/, $_);
		$family = $a[1];
		
		if (defined $itak_obj{$cotyledon}{$species}{$family}) {
			$itak_obj{$cotyledon}{$species}{$family}.="\t".$a[0];
		} else {
			$itak_obj{$cotyledon}{$species}{$family} = $a[0];
		}

		if (defined $protein_family{$a[0]}) {
			$protein_family{$a[0]}.="\t".$a[1];
		} else {
			$protein_family{$a[0]} = $a[1];
		}

		$protein_species{$a[0]} = $species;
	}
	$pcat->close;
	$tcat->close;

	# transcripts id and gene id to hash
	# protein id and transcript id to hash
	if (-s $trans_gene_file)
	{
		my $tgf = IO::File->new($trans_gene_file) || die "Can not open trans gene file $!\n";
		while(<$tgf>)
		{
			chomp;
			my ($protein_id, $transcript_id, $gene_id) = split(/\t/, $_);
			if (defined $trans_gene{$transcript_id}) 
			{ 
				die "Error, transcrips $transcript_id has two corresponding gene IDs :".
				    " $trans_gene{$transcript_id} and $gene_id\n"; 
			}
			else 
			{
				$trans_gene{$transcript_id} = $gene_id;
			}

			if ( defined $trans_protein{$transcript_id} )
			{
				die "Error, transcript id $transcript_id has two protein ID : $trans_protein{$transcript_id} $protein_id \n";
			}
			else 
			{
				$trans_protein{$transcript_id} = $protein_id;
			}

			if ( defined $protein_trans{$protein_id} )
			{
				die "Error, protein id $protein_id has two transcript ID : $protein_trans{$protein_id} $transcript_id \n";
			}
			else
			{
				$protein_trans{$protein_id} = $transcript_id;
			}
		}
		$tgf->close;
	}

	# species category to hash;
	$species_category{$species} = $cotyledon;

	# sequence to hash
	my $pin = Bio::SeqIO->new(-format=>'fasta', -file=>$pk_seq);
	while(my $pkseq = $pin->next_seq)
	{
		$protein_seq{$pkseq->id} = $pkseq->seq;
	}

	my $tin = Bio::SeqIO->new(-format=>'fasta', -file=>$tf_seq);
	while(my $tfseq = $tin->next_seq)
	{
		$protein_seq{$tfseq->id} = $tfseq->seq;
	}

	# gene domain to hash
	my $gd1 = IO::File->new($pk_aln) || die "Can not open PKs alignment file $pk_aln $!\n";
	while(<$gd1>)
	{
		chomp;
		my @a = split(/\t/, $_);	
		$protein_domain{$a[0]."\t".$a[1]} = $a[2]."\t".$a[3]."\t".$a[4]."\t".$a[5]."\t".$a[6]."\t".$a[7]."\t".$a[8]."\t".$a[9]."\t".$a[10]; 
	}
	$gd1->close;

	my $gd2 = IO::File->new($tf_aln) || die "Can not open TFs alignment file $tf_aln $!\n";
	while(<$gd2>)
	{
		chomp;
		my @a = split(/\t/, $_);
		$protein_domain{$a[0]."\t".$a[1]} = $a[2]."\t".$a[3]."\t".$a[4]."\t".$a[5]."\t".$a[6]."\t".$a[7]."\t".$a[8]."\t".$a[9]."\t".$a[10];
	}
	$gd2->close;

	# cds sequence to hash
	my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$cds_file);
	while(my $inseq = $in->next_seq)
	{
		my $trans_id = $inseq->id;

		unless (defined $trans_protein{$trans_id} ) { die "Error in transcript id $trans_id, do not have protein ID $species\n"; }

		my $protein_id = $trans_protein{$trans_id};
		
		if ( defined $protein_seq{$protein_id} )
		{
			$transcript_seq{$inseq->id} = $inseq->seq;
		}
	}
}
$fh->close;

#================================================================
# prepare files for cluster / trees
#================================================================

my %all_specie;
my %all_family;

foreach my $cotyledon (sort keys %itak_obj)
{
	foreach my $species (sort keys %{$itak_obj{$cotyledon}})
	{

		$all_specie{$species} = 1;

		foreach my $family (sort keys %{$itak_obj{$cotyledon}{$species}})
		{
			# get pep sequence base on cluster family
			my ($protein_seq1, $protein_seq2) = ("", "", "");
			
			my $main_gene_id; my %main_gene;
			my @proteins = split(/\t/, $itak_obj{$cotyledon}{$species}{$family});
			foreach my $protein_id (@proteins) {
				unless (defined $protein_trans{$protein_id}) { die "Error, do not have transcript id for $protein_id\n"; }
				my $transcript_id = $protein_trans{$protein_id};
				$main_gene_id = $trans_gene{$transcript_id};
				$main_gene{$main_gene_id} = 1;
			}

			my $type;
			if ($family =~ m/^PPC/) { $type = "pk"; } else { $type = "tf"; }


			$all_family{$family} = 1;

			$family_sum{$species}{$family} = scalar(keys(%main_gene));
		}
	}
}

print "Name";
foreach my $fam_name (sort keys %all_family) { print "\t".$fam_name; }
print "\n";

foreach my $species (sort keys %all_specie)
{
	print $species;

        foreach my $fam_name (sort keys %all_family)
        {
		if ( defined $family_sum{$species}{$fam_name} )
		{
                	print "\t".$family_sum{$species}{$fam_name};
		}
		else
		{
			print "\t0";
		}
        }
	print "\n";
}

