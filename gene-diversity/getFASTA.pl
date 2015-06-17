#!/usr/bin/perl

#--------------------------------------------------------
# PROGRAM: getFASTA.pl
# AUTHOR:  Sofia Hauck
# CREATED: 09.06.2015
# UPDATED: ----------
# VERSION: v1.00
#--------------------------------------------------------
# VERSION HISTORY
# v1.00 (09.06.2015) created
#--------------------------------------------------------

use strict;
use warnings;
use LWP::Simple;

# Declares subroutines
sub Usage( ; $ );

# Get Command line options, exits if conditions don't look right
if( scalar(@ARGV) < 1 ) { Usage("Not enough command line options"); exit; }
my $fIn = $ARGV[0];

# Command line option checks
if(! defined $fIn) { Usage("Missing Option: input file <FILE>"); exit; }
if(! -e $fIn) { Usage("Input file does not exist: $fIn"); exit; }

# Preparing for output
my $dOut = $fIn . "-files";
if(-e $dOut) { Usage("Output directory already exist: $dOut"); exit; }
mkdir $dOut; 

# Open infile, move to outfile directory
open(INFILE, $fIn) or die "Cannot open $fIn\n";
	chdir $dOut;
	my $database = <INFILE>;
	chomp $database;
	
	# For each locus in the list, get URL of FASTA file and save it to directory
	while ( my $line = <INFILE> )											
	{
		chomp $line; 
		
		my $url = "http://rest.pubmlst.org/db/".$database."/loci/".$line."/alleles_fasta";
		print "URl is $url \n";
		
		my $file = $line.".FAS";
		
		getstore($url, $file);
	}
close(INFILE);

print "Finished! Results are in $dOut directory.\n";

#---------------------------------------------------------------
# Subroutines
#---------------------------------------------------------------

sub Usage( ; $ )
{
	my $message = $_[0] || '';

	print << 'EOU';

getFASTA.pl

Description:
  From a text file with a list of loci (one per column), gets all the FASTA files via RESTful API. 
  Currently set up to work only with mycobacteria_seqdef database. 
	  
Usage:
getFASTA.pl [ options ]

Results will be in folder with same name as file,
with files within the folder being named after the locus.

EOU

print "Quit because: $message\n";
print "ARGV was " . join (", ", @ARGV) . "\n";
}

