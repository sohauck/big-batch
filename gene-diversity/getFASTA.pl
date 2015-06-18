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

$| =1; # for dynamic output

# Declares subroutines
sub Usage( ; $ );

# Get Command line options, exits if conditions don't look right
my $fIn;  # file with the database name and loci names
my $dOut; # directory where exported FASTA files will go

if( scalar(@ARGV) < 2 ) { Usage("Not enough command line options"); exit; }
my $i = 0; 
my $arg_cnt = 0; 
for ($i=0; $i<=$#ARGV; $i++)
{
	if($ARGV[$i] eq "-h")	       { Usage("You asked for help"); exit; }
	if($ARGV[$i] eq "-in")         { $fIn  = $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-out")        { $dOut = $ARGV[$i+1] || ''; $arg_cnt++; }
}

# Command line option checks
if(! defined $fIn) { Usage("Missing Option: input file <FILE>"); exit; }
if(! -e $fIn) { Usage("Input file does not exist: $fIn"); exit; }

# Preparing for folder for output
if( -e $dOut) { Usage("Output directory already exist: $dOut"); exit; }
mkdir $dOut; 

print "Extraction now up to...\n";

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
		
		my $file = $line.".FAS";
		getstore($url, $file);
		
		print "\r$line";
	}
close(INFILE);

print "\nFinished! Results are in $dOut directory.\n";

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
 -in FILE with database in first line, loci in the rest
 -out DIRECTORY file path including new folder name, for exported files to go into



EOU

print "Quit because: $message\n";
print "ARGV was " . join (", ", @ARGV) . "\n";
}

