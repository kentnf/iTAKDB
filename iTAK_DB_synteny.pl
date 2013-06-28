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
  1) file prefix
  2) cotyledon ( monocotyledon, dicotyledon, non-angiosperms )
  3) species

  * input files base on file prefix
  1) proteins
  2) CDS
  3) transcript gene file
  4) iTAK output

  * the prefix is for below files
  0) prefix_pep (get prefix from input protein sequence)
  1) prefix_chrSize
  2) prefix_gene_position
  3) prefix_block.csv
  4) prefix_trans_gene

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
my %species;

my $cfh = IO::File->new($list) || die "Can not open list file $list $!\n";
while(<$cfh>)
{
	chomp; 
	if ($_ =~ m/^#/) { next; }

	my ($ipath, $cotyledon, $species, 
	    $pep_file, $trans_gene_file, $pk_cat, $tf_cat,
	    $block, $chrSize, $gene_position,
	    $ref_prefix, $ref_block);

	my @a = split(/\t/, $_);
	if ( scalar(@a) == 3 )		{ ($ipath, $cotyledon, $species) = @a; }
	elsif ( scalar(@a) == 5 ) 	{ ($ipath, $cotyledon, $species, $ref_prefix, $ref_block) = @a; } 
	else 				{  die "Error in list file $_ \n"; } 

	# species to hash and check the repeat species
	if (defined $species{$species}) { die "Error, repeat species $species\n"; } else { $species{$species} = 1;}

	# set files
	my ($prefix,$rpath,$rsuffix) = fileparse($ipath);
	$pep_file 	= $ipath."/".$prefix."_pep";
	$trans_gene_file= $ipath."/".$prefix."_trans_gene";
	$chrSize 	= $ipath."/".$prefix."_chrSize";
	$gene_position 	= $ipath."/".$prefix."_gene_position";
	$block 		= $ipath."/".$prefix."_block";

	# set iTAK output files
	# output folder is: TAIR9_protein_output
	# output files are: TAIR9_protein_pkaln, TAIR9_protein_pkcat, TAIR9_protein_pkseq
	#                   TAIR9_protein_tf_align, TAIR9_protein_tf_family, TAIR9_protein_tf_seq
	my $folder = $pep_file."_output";
	my ($fname,$fpath,$fsuffix) = fileparse($pep_file);
	$pk_cat = $folder."/".$fname."_pkcat";
	$tf_cat = $folder."/".$fname."_tf_family";

	# check the input files
	my @files = ($pk_cat, $tf_cat, $chrSize, $gene_position, $block, $trans_gene_file);
	foreach my $file ( @files ) {
		unless (-s $file) { die "Error! File $file do not exist or has content $!\n"; }
	}

	if ($ref_block) {
		unless (-s $ref_block) { die "Error! File reference block $ref_block do not exist $!\n"; }
	}
}
$cfh->close;

#================================================================
# list file into to hash					
#================================================================
# init hash for clustering : key1: cyto; key2: species; key3: TF or PK; value: genes (tab delimit)
my %itak_obj;

# hash for synteny file
my %prefix_species;
my %prefix_id_type;
my %prefix_chrSize;
my %prefix_gene_position;
my %prefix_block;
my %prefix_ref;
my %prefix_ref_block;

# has for protein / transcript / gene ID
my %trans_gene;		# key: transcript_id, value: gene_id
my %trans_protein;	# key: transcript_id, value: protein_id
my %protein_trans;	# key: protein_id, value: transcript_id

my $fh = IO::File->new($list) || die "Can not open list file $list $!\n";
while(<$fh>)
{
	chomp;
	if ($_ =~ m/^#/) { next; }

	my ($ipath, $cotyledon, $species,
	    $pep_file, $trans_gene_file, $pk_cat, $tf_cat,
	    $block, $chrSize, $gene_position,
	    $ref_prefix, $ref_block);

	my @a = split(/\t/, $_);
	if ( scalar(@a) == 3 )          { ($ipath, $cotyledon, $species) = @a; }
	elsif ( scalar(@a) == 5 )       { ($ipath, $cotyledon, $species, $ref_prefix, $ref_block) = @a; }
	else                            {  die "Error in list file $_ \n"; }

	# species to hash and check the repeat species
	if (defined $species{$species}) { die "Error, repeat species $species\n"; } else { $species{$species} = 1;}

	# set files
	my ($prefix,$rpath,$rsuffix) = fileparse($ipath);
	$pep_file       = $ipath."/".$prefix."_pep";
	$trans_gene_file= $ipath."/".$prefix."_trans_gene";
	$chrSize        = $ipath."/".$prefix."_chrSize";
	$gene_position  = $ipath."/".$prefix."_gene_position";
	$block          = $ipath."/".$prefix."_block";

	# set iTAK output file
	# output folder is: TAIR9_protein_output
	# output files are: TAIR9_protein_pkaln, TAIR9_protein_pkcat, TAIR9_protein_pkseq
	#                   TAIR9_protein_tf_align, TAIR9_protein_tf_family, TAIR9_protein_tf_seq
	my $folder = $pep_file."_output";
        my ($fname,$fpath,$fsuffix) = fileparse($pep_file);
        $pk_cat = $folder."/".$fname."_pkcat";
        $tf_cat = $folder."/".$fname."_tf_family";

	$prefix_species{$prefix} = $species;
	$prefix_chrSize{$prefix} = $chrSize;
	$prefix_gene_position{$prefix} = $gene_position;
	$prefix_block{$prefix} = $block;

	if ($ref_prefix) 
	{ 
		$prefix_ref{$prefix} = $ref_prefix; 
		$prefix_ref_block{$prefix} = $ref_block;
	}

	# classification to hash
	my $family;
	my $pcat = IO::File->new($pk_cat) || die "Can not open pk cat file $pk_cat $!\n";
	my $tcat = IO::File->new($tf_cat) || die "Can not open tf cat file $tf_cat $!\n";
	while(<$pcat>) 
	{
		chomp;
		my @a = split(/\t/, $_);
		$family = $a[1];

		if (defined $itak_obj{$cotyledon}{$prefix}{$family}) {
			$itak_obj{$cotyledon}{$prefix}{$family}.="\t".$a[0];
		} else {
			$itak_obj{$cotyledon}{$prefix}{$family} = $a[0];
		}
	}
	while(<$tcat>)
	{
		chomp;
		my @a = split(/\t/, $_);
		$family = $a[1];
		
		if (defined $itak_obj{$cotyledon}{$prefix}{$family}) {
			$itak_obj{$cotyledon}{$prefix}{$family}.="\t".$a[0];
		} else {
			$itak_obj{$cotyledon}{$prefix}{$family} = $a[0];
		}
	}
	$pcat->close;
	$tcat->close;

	# transcripts id and gene id to hash
	# protein id and transcript id to hash
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
$fh->close;

#================================================================
# prepare files for synteny
#================================================================
mkdir("for_synteny");

foreach my $cotyledon (sort keys %itak_obj)
{
	foreach my $prefix (sort keys %{$itak_obj{$cotyledon}})
	{
		# get files for plant synteny 
		my ($chrSize, $gene_position, $block);
		$chrSize = $prefix_chrSize{$prefix};
		$gene_position = $prefix_gene_position{$prefix}; 
		$block = $prefix_block{$prefix};
		
		# get files for plant and reference synteny
		my ($ref_prefix, $ref_chrSize, $ref_gene_position, $ref_block);
		if ( defined $prefix_ref{$prefix} ) 
		{ 
			$ref_prefix = $prefix_ref{$prefix}; 
			$ref_chrSize = $prefix_chrSize{$ref_prefix};
			$ref_gene_position = $prefix_gene_position{$ref_prefix};
			$ref_block = $prefix_ref_block{$prefix};
		}

		# get species name
		unless (defined $prefix_species{$prefix}) { die "Error, do not find species name for $prefix $!\n"; }
		my $species = $prefix_species{$prefix};

		# parse gene/mRNA ID for each family, then generate synteny
		foreach my $family (sort keys %{$itak_obj{$cotyledon}{$prefix}})
		{	
			my $gene_id; my %uniq_gene; my $trans_id; my %uniq_trans; my %uniq_protein;
			my @proteins = split(/\t/, $itak_obj{$cotyledon}{$prefix}{$family});
			foreach my $protein_id (@proteins) 
			{	
				$uniq_protein{$protein_id} = 1;
	
				$trans_id = $protein_trans{$protein_id};
				$uniq_trans{$trans_id} = 1;

				$gene_id = $trans_gene{$trans_id};
				$uniq_gene{$gene_id} = 1;
			}

			# save protein family info to file
			my $type;
			if ($family =~ m/^PPC/) { $type = "pk"; } else { $type = "tf"; }
			my $fm = $family;  $fm =~ s/\//_/;  $fm =~ s/\s+/_/;  $fm =~ s/:/_/;
			my $sp = $species; $sp =~ s/\s+/_/;

                        my $syn_cmd = "";
                        my $gene_family_list = "for_synteny/".$sp."_".$fm;
                        my $fh3 = IO::File->new(">".$gene_family_list) || die "Can not open gene family list file $gene_family_list $!\n";
                        foreach my $id (sort keys %uniq_gene) { print $fh3 "$id\t$id\n"; }
                        $fh3->close;

                        if ( defined $prefix_ref{$prefix} ) {
                        	$syn_cmd = "perl ../synteny/plant_synteny.pl -i $gene_family_list -a $gene_position -b $block -c $chrSize -x $ref_gene_position -y $ref_block -z $ref_chrSize";
                                print $syn_cmd."\n";
                        }
			else
			{
				$syn_cmd = "perl ../synteny/plant_synteny.pl -i $gene_family_list -a $gene_position -b $block -c $chrSize";
				print $syn_cmd."\n";
			}
		}
	}
}

#================================================================
# kentnf: subroutine						
#================================================================
