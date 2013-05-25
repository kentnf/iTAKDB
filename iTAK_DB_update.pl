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

#===============================================================#
# PKs ID and Description to hash				#
#===============================================================#
my $pk_desc_file = ${FindBin::RealBin}."/protein_kinase_family";
unless (-s $pk_desc_file) { die "Can not locate protein_kinase_family file\n"; }
my %pk_desc;	# key: PKsID, value:desc;
my $kfh = IO::File->new($pk_desc_file) || die "Can not open file $pk_desc_file $!\n";
while(<$kfh>)
{
	chomp;
	my @a = split(/\t/, $_);
	$pk_desc{$a[0]} = $a[1];
}
$kfh->close;

#================================================================
# prepare files for cluster / trees
#================================================================
mkdir("for_cluster");
mkdir("for_tree");

foreach my $cotyledon (sort keys %itak_obj)
{
	foreach my $species (sort keys %{$itak_obj{$cotyledon}})
	{
		foreach my $family (sort keys %{$itak_obj{$cotyledon}{$species}})
		{
			# get pep sequence base on cluster family
			my ($protein_seq1, $protein_seq2) = ("", "", "");
			
			my $main_gene_id; my %main_gene;
			my @proteins = split(/\t/, $itak_obj{$cotyledon}{$species}{$family});
			foreach my $protein_id (@proteins) {
				unless (defined $protein_seq{$protein_id}) { die "Error, do not have seq for $protein_id\n"; }
				$protein_seq1 .= ">".$protein_id."\n".$protein_seq{$protein_id}."\n";
				$protein_seq2 .= ">".$protein_id."\n".$protein_seq{$protein_id}."\n";
				
				unless (defined $protein_trans{$protein_id}) { die "Error, do not have transcript id for $protein_id\n"; }
				my $transcript_id = $protein_trans{$protein_id};
				$main_gene_id = $trans_gene{$transcript_id};
				$main_gene{$main_gene_id} = 1;
			}

			if ($cotyledon eq "dicotyledon") {
				if ($species ne "Arabidopsis") {
					if ( defined $itak_obj{"dicotyledon"}{"Arabidopsis"}{$family} )
					{
						my @di_proteins = split(/\t/, $itak_obj{"dicotyledon"}{"Arabidopsis"}{$family});
						foreach my $protein_id ( @di_proteins ) {
							unless (defined $protein_seq{$protein_id}) { die "Error, do not have seq for $protein_id\n"; }
							$protein_seq1 .= ">".$protein_id."\n".$protein_seq{$protein_id}."\n";
						}
					}
				}
			} elsif ( $cotyledon eq "monocotyledon") {
				if ($species ne "Rice") {
					if ( defined $itak_obj{"monocotyledon"}{"Rice"}{$family} )
					{
						my @mono_proteins = split(/\t/, $itak_obj{"monocotyledon"}{"Rice"}{$family});
						foreach my $protein_id (@mono_proteins) {
							unless (defined $protein_seq{$protein_id}) { die "Error, do not have seq for $protein_id\n"; }
							$protein_seq1 .= ">".$protein_id."\n".$protein_seq{$protein_id}."\n";
						}
					}
				}
			} elsif ( $cotyledon eq "non-angiosperms") {
				# code for generate pep sequence for each family of non-angiosperms
			} else {
				die "Error at cotyledon $cotyledon $!\n";
			}

			# save pep sequences to file
			my $type;
			if ($family =~ m/^PPC/) { $type = "pk"; } else { $type = "tf"; }
			my $fm = $family; 
			$fm =~ s/\//_/;  
			$fm =~ s/\s+/_/;
			$fm =~ s/:/_/;
			my $sp = $species;
			$sp =~ s/\s+/_/;

			#my $cluster_file = "for_cluster/".$cotyledon."/".$type."/".$sp."_".$fm.".pep";
			my $cluster_file = "for_cluster/".$sp."_".$fm.".pep";
			my $tree_file    = "for_tree/".$sp."_".$fm.".pep";

			my $fh1 = IO::File->new(">".$cluster_file) || die "Can not open cluster file $cluster_file $!\n";
			print $fh1 $protein_seq1;
			$fh1->close;

			my $fh2 = IO::File->new(">".$tree_file) || die "Can not open tree file $tree_file $!\n";
			print $fh2 $protein_seq2;
			$fh2->close;
			
			cluster($cluster_file);
			cluster($tree_file);

			# number of gene in each family
			$family_sum{$species}{$family} = scalar(keys(%main_gene));
			#print $cotyledon,"\t",$species,"\t",$family,"\t",scalar(keys(@genes))."\n";
		}
	}
}

#===============================================================#
# generate files for download					#
#===============================================================#
mkdir("for_download");
my ($all_TFs_pep, $all_TFs_cds, $all_PKs_pep, $all_PKs_cds, $all_TFs_list, $all_PKs_list) = (
	"for_download/All_TFs_TRs_protein", 
	"for_download/All_TFs_TRs_nucleotide",
	"for_download/All_PKs_protein",
	"for_download/All_PKs_nucleotide",
	"for_download/All_TFs_TRs_list",
	"for_download/All_PKs_list"
);

my ($dfh1, $dfh2, $dfh3, $dfh4, $dfh5, $dfh6);
$dfh1 = IO::File->new(">".$all_TFs_pep) || die "Can not open file $all_TFs_pep $!\n";
$dfh2 = IO::File->new(">".$all_TFs_cds) || die "Can not open file $all_TFs_cds $!\n";
$dfh3 = IO::File->new(">".$all_PKs_pep) || die "Can not open file $all_PKs_pep $!\n";
$dfh4 = IO::File->new(">".$all_PKs_cds) || die "Can not open file $all_PKs_cds $!\n";
$dfh5 = IO::File->new(">".$all_TFs_list)|| die "Can not open file $all_TFs_list $!\n";
$dfh6 = IO::File->new(">".$all_PKs_list)|| die "Can not open file $all_PKs_list $!\n";

my %protein_cds_ID;

foreach my $protein_id (sort keys %protein_family)
{
	my $transcript_id;
 	if (defined $protein_trans{$protein_id} )
	{
		$transcript_id = $protein_trans{$protein_id};
	}
	else
	{
		die "Error, protein $protein_id do not have corresponding transcript id\n";
	}

	my @family = split(/\t/, $protein_family{$protein_id});
	foreach my $family (@family) 
	{
		if ($family =~ m/^PPC/) 
		{
			print $dfh3 ">".$protein_id."\t".$family." ".$pk_desc{$family}."\n".$protein_seq{$protein_id}."\n";
			print $dfh4 ">".$transcript_id."\t".$family." ".$pk_desc{$family}."\n".$transcript_seq{$transcript_id}."\n";
			print $dfh6 $protein_id."\t".$family."\t".$pk_desc{$family}."\n";
		}
		else 
		{
			print $dfh1 ">".$protein_id."\t".$family."\n".$protein_seq{$protein_id}."\n";
			print $dfh2 ">".$transcript_id."\t".$family."\n".$transcript_seq{$transcript_id}."\n";
			print $dfh5 $protein_id."\t".$family."\n";
		}
	}
}

$dfh1->close;
$dfh2->close;
$dfh3->close;
$dfh4->close;
$dfh5->close;
$dfh6->close;

my ($tf_folder, $pk_folder) = ("for_download/transcription_factors_and_transcriptional_regulators", "for_download/protein_kinases");
mkdir($tf_folder); mkdir($pk_folder);

my ($species_PKs_pep, $species_PKs_cds, $species_TFs_pep, $species_TFs_cds);

foreach my $cotyledon (sort keys %itak_obj)
{
        foreach my $species (sort keys %{$itak_obj{$cotyledon}})
        {
		print $species."\n";
		my $sp = $species;
		$sp =~ s/\s+/_/;
		($species_PKs_pep, $species_PKs_cds, $species_TFs_pep, $species_TFs_cds) = (
		$pk_folder."/".$sp."_protein",
		$pk_folder."/".$sp."_cds",
		$tf_folder."/".$sp."_protein",
		$tf_folder."/".$sp."_cds");

		my $pk_pep_fh = IO::File->new(">".$species_PKs_pep) || die "Can not open file pk pep for $species $!\n";
		my $pk_cds_fh = IO::File->new(">".$species_PKs_cds) || die "Can not open file pk cds for $species $!\n";
		my $tf_pep_fh = IO::File->new(">".$species_TFs_pep) || die "Can not open file tf pep for $species $!\n";
		my $tf_cds_fh = IO::File->new(">".$species_TFs_cds) || die "Can not open file tf cds for $species $!\n";

		my ($tf_pep_seq, $tf_cds_seq, $pk_pep_seq, $pk_cds_seq);

                foreach my $family (sort keys %{$itak_obj{$cotyledon}{$species}})
                {
			my @proteins = split(/\t/, $itak_obj{$cotyledon}{$species}{$family});
			foreach my $protein_id ( @proteins ) {
				my $transcript_id = $protein_trans{$protein_id};
				unless(defined $protein_seq{$protein_id}) { die "Error, do not find protein of $protein_id in $species \n"; }
				unless(defined $transcript_seq{$transcript_id}) { die "Error, do not find CDS of $transcript_id in $species \n"; }
				if ($family =~ m/^PPC/) {
					$pk_pep_seq.=">".$protein_id."\t".$family." ".$pk_desc{$family}."\n".$protein_seq{$protein_id}."\n";
					$pk_cds_seq.=">".$transcript_id."\t".$family." ".$pk_desc{$family}."\n".$transcript_seq{$transcript_id}."\n";
				} else {
			

					$tf_pep_seq.=">".$protein_id."\t".$family."\n".$protein_seq{$protein_id}."\n";
					$tf_cds_seq.=">".$transcript_id."\t".$family."\n".$transcript_seq{$transcript_id}."\n";
				}
			}
		}

		print $pk_pep_fh $pk_pep_seq;
	        print $pk_cds_fh $pk_cds_seq;
	        print $tf_pep_fh $tf_pep_seq;
	        print $tf_cds_fh $tf_cds_seq;

		$pk_pep_fh->close;
		$pk_cds_fh->close;
		$tf_pep_fh->close;
		$tf_cds_fh->close;
	}
}

#================================================================
# generate files for blast					
#================================================================
mkdir("for_blast");
my ($blast_TFs_cds, $blast_TFs_pep, $blast_PKs_cds, $blast_PKs_pep) = (
        "for_blast/All_TFs_TRs_nucleotide",
        "for_blast/All_TFs_TRs_protein",
        "for_blast/All_PKs_nucleotide",
        "for_blast/All_PKs_protein"
);

parse_seq_for_blast($all_TFs_pep, $blast_TFs_pep);
parse_seq_for_blast($all_TFs_cds, $blast_TFs_cds);
parse_seq_for_blast($all_PKs_pep, $blast_PKs_pep);
parse_seq_for_blast($all_PKs_cds, $blast_PKs_cds);
#system("formatdb -i $blast_TFs_pep -p T");
#system("formatdb -i $blast_TFs_cds -p F");
#system("formatdb -i $blast_PKs_pep -p T");
#system("formatdb -i $blast_PKs_cds -p F");

#================================================================
# generate files and cmds for load data to database.
#================================================================

# create folder
my $db_folder = "for_database";
mkdir($db_folder);

# step 1  generate family_table file
# do not need generate this table unless new family was add to itak
#my $fh1 = IO::File->new(">family_table") || die "Can not open file family_table $!\n";
#foreach my $fam_name (sort keys %family)
#{
#	print $fam_name."\t".$family{$fam_name}."\n";
#}
#$fh1->closel

#================================================================
# step2 generate family_summary file
# family	species		number of genes
# ABI3VP1	arabidopsis	72
#================================================================
my $family_summary = $db_folder."/family_summary";
my $fh2 = IO::File->new(">".$family_summary) || die "Can not open file family_summary $!\n";
foreach my $species (sort keys %family_sum)
{
	foreach my $fam_name (sort keys %{$family_sum{$species}})
	{
		print $fh2 $fam_name."\t".$family_sum{$species}{$fam_name}."\t".$species."\n";
	}
}
$fh2->close;

#================================================================
# step3 generate gene_family file
# gene		family
# AT1G01010.1	NAC
#================================================================
my $protein_family_table = $db_folder."/protein_family_table";
my $fh3 = IO::File->new(">".$protein_family_table) || die "Can not open file protein_family_table $!\n";
foreach my $protein_id (sort keys %protein_family)
{
	my @family = split(/\t/, $protein_family{$protein_id});
	foreach my $family (@family)
	{
		print $fh3 $protein_id."\t".$family."\n";
	}
}
$fh3->close;

#================================================================
# step4 generate gene_domain file
#================================================================
my $protein_domain_table = $db_folder."/protein_domain_table";
my $fh4 = IO::File->new(">".$protein_domain_table) || die "Can not open file protein_domain_talbe $!\n";
foreach my $key (sort keys %protein_domain)
{
	print $fh4 $key."\t".$protein_domain{$key}."\n";
}
$fh4->close;

#================================================================
# step5 generate gene_table file
#================================================================
my $gene_table = $db_folder."/gene_table";
my $fh5 = IO::File->new(">".$gene_table) || die "Can not open file gene_table $!\n";
foreach my $protein_id (sort keys %protein_seq)
{
	unless ( defined $protein_trans{$protein_id} ) { die;  }
	my $transcript_id = $protein_trans{$protein_id};
	unless ( defined $trans_gene{$transcript_id}) { die; }
	my $gene_id = $trans_gene{$transcript_id};

	print $fh5 $protein_id."\t".$gene_id."\t".$protein_seq{$protein_id}."\t".$transcript_seq{$transcript_id}."\t".
		   $protein_species{$protein_id}."\n";
}
$fh5->close;

#================================================================
# step6 generate category table					
#================================================================
my $category_table = $db_folder."/category";
my $fh6 = IO::File->new(">".$category_table) || die "Can not open file category $!\n";
foreach my $species (sort keys %species_category)
{
	print $fh6 $species."\t".$species_category{$species}."\n";
}
$fh6->close;

# some info about update database
print "Load data to mysql database......\n";
#print "load data local infile \"domain_table\" into table domain;\n";
#print "load data local infile \"family_table\" into table family;\n";
print "load data local infile \"family_summary\" into table family_summary;\n";
print "load data local infile \"gene_table\" into table gene;\n";
print "load data local infile \"gene_annotation\" into table gene_annotation;\n";
print "load data local infile \"protein_domain_table\" into table gene_domain;\n";
print "load data local infile \"protein_family_table\" into table gene_family;\n";

#================================================================
# kentnf: subroutine						
#================================================================
=head1

 parse_seq_for_blast

=cut
sub parse_seq_for_blast
{
	my ($in_file, $out_file) = @_;

	my $out = IO::File->new(">".$out_file) || die "Can not open output file $out_file $!\n";
	my $in = IO::File->new($in_file) || die "Can not open input file $in_file $!\n";
	while(<$in>)
	{
		$_ =~ s/\t/ /ig;
		print $out $_;
	}
	$in->close;
	$out->close;
}

=head1 

 cluster

=cut
sub cluster
{
	my $protein_seq = shift;

	my $seq_num;
	my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$protein_seq);
	while(my $inseq = $in->next_seq) { $seq_num++; }

	my $nwk_file = $protein_seq; $nwk_file =~ s/\.pep/\.nwk/;
	my $aln_file = $protein_seq; $aln_file =~ s/\.pep/\.aln/;
	my $phb_file = $protein_seq;

	if ( $seq_num >= 3 )
	{
		if ( $seq_num > 3 )
		{
			unless ( -s $nwk_file )
			{
				$phb_file =~ s/\.pep/\.phb/;
				my $cmd1 = "$clustalw2_program -INFILE=$protein_seq -ALIGN";
				my $cmd2 = "$clustalw2_program -INFILE=$protein_seq -PROFILE1=$aln_file -BOOTSTRAP=100";
				system($cmd1) && die "Error in command $cmd1\n";
                                system($cmd2) && die "Error in command $cmd2\n";
                                unless (-s $phb_file) { die "Error, do not have phb file $phb_file\n"; }
                                phb2nwk($phb_file, $nwk_file);
			}
		}
		else
		{
			unless ( -s $nwk_file )
			{
				$phb_file =~ s/\.pep/\.ph/;
				my $cmd1 = "$clustalw2_program -INFILE=$protein_seq -TREE";
				system($cmd1) && die "Error in command $cmd1\n";
				unless (-s $phb_file) { die "Error, do not have phb file $phb_file\n"; }
				phb2nwk($phb_file, $nwk_file);
			}
		}
	}
}

=head1

 phb2nwk -- convert phb/ph format to nwk format

=cut
sub phb2nwk
{
        my ($phb, $nwk) = @_;

        my $fh_phb = IO::File->new($phb) || die "Can not open phb file: $phb \n";
        my $fh_nwk = IO::File->new(">".$nwk) || die "Can not open nwk file: $nwk \n";

        while(<$fh_phb>)
        {
                chomp;
                if ($_ =~ m/^:(\d+)\.(\d+)\[(\d+)\](.*)/)
                {
                        print $fh_nwk $3.":".$1.".".$2.$4;
                }
                elsif ($_ =~ m/^(\S+):(\d+)\.(\d+)(.*)/)
                {
                        print $fh_nwk $_;
                }
                elsif ($_ =~ m/^(\S+):-(\d+)\.(\d+)(.*)/)
                {
                        print $fh_nwk $1.":0.000".$4;
                }
                else
                {
                        print $fh_nwk $_;
                }
        }

        print $fh_nwk $_;

        $fh_phb->close;
        $fh_nwk->close;
}
