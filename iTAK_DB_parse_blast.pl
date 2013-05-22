#!/usr/bin/perl

=head

 parse blast raw result to tab delimit table

=cut

use strict;
use warnings;
use Bio::SearchIO;
use Getopt::Long;

my $usage = "\nusage: $0 blast_result > output\n\n";
my $input_file = shift || die $usage;

my $max_hit_num ||= 5;

#===============================================================#
# main								#
#===============================================================#
my ($hit_num, $query_name, $query_length);

my $in = Bio::SearchIO->new(-format=>'blast', -file=>$input_file);
while(my $result = $in->next_result)
{
	$hit_num = 0;

	$query_name = $result->query_name;
	$query_length = $result->query_length;

	while(my $hit = $result->next_hit)
	{
		if ($query_name ne $hit->name)
		{
			$hit_num++;

			if ($hit_num <= $max_hit_num)
			{
				print 	$query_name."\t".
				      	$hit->accession."\t".
					$hit->description."\t".
					$hit->raw_score."\t".
					$hit->significance."\n";
			}
		}
	}
}
