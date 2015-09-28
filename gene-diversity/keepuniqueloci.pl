#!/usr/bin/perl

#--------------------------------------------------------
# PROGRAM: keepuniqueloci.pl
# AUTHOR:  Sofia Hauck
# CREATED: 16.09.2015
#--------------------------------------------------------

use strict;
use warnings;
$| = 1;

# Declares subroutines
sub Usage( ; $ );

# Defines scalars needed from command line
my $fTab;    # file with isolate/locus table where relevant loci are named
my $dIn;    # directory where all the FASTA files for each locus are
my $dOut;    # directory where only the relevant parts of those FASTA files will be copied to
my $dup = 1; #  turned to the smallest number of duplicates necessary for locus to be considered, default is 1

# Get Command line options, exits if conditions don't look right
if( scalar(@ARGV) < 1 ) { Usage("Not enough command line options"); exit; }
my $i = 0;
my $arg_cnt = 0; 
for ($i=0; $i<=$#ARGV; $i++)
{
	if($ARGV[$i] eq "-h")	          { Usage("You asked for help"); exit; }
	if($ARGV[$i] eq "-din")           { $dIn = $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-dout")          { $dOut = $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-dup")           { $dup = $ARGV[$i+1] || ''; $arg_cnt++; }
}

# Command line option checks & File checks
if(! defined $dIn)  { Usage("Missing Option: -din DIRECTORY"); exit; }
if(! defined $dOut)  { Usage("Missing Option: -dout DIRECTORY"); exit; }

if(! -e $dIn)  { Usage("Input directory doesn't exist: $dIn"); exit; }
if(  -e $dOut)  { Usage("Output directory already exists: $dOut"); exit; }

# Makes output directory and prepared output count file, including header
mkdir $dOut;

my $uniquefile = $dOut . "-count-nuc.txt";

if( -e $uniquefile)  { Usage("Output file already exists: $uniquefile"); exit; }

open(UNIQUENUC, '>', $uniquefile) or die "Cannot open $uniquefile\n";
	print UNIQUENUC "locus,count-nuc\n";


# Get names of all the files in the directory
opendir (ORIGDIR, $dIn) or die "Cannot open directory: $!";
	my @files = readdir ORIGDIR;
	@files = grep(/^([A-Z]|[a-z]|[0-9])/,@files);
	if ($#files < 1)
	{ die "It looks like you have no FASTA files to align. If their file names don't begin with letters or numbers they are not being included.\n";}
closedir ORIGDIR;


# Loop that removes duplicates and renumbers alleles in all the files and puts them in new directory

foreach my $file (@files)
{
 	my $infile =  $dIn  . "/" . $file; # original file from folder given in command line
 	my $outfile = $dOut . "/" . $file; # aligned file into translated folder, with same file name
 	
	my %unique = (); # empty hash to for uniqueness 
	my $allelenumber = 1;
	my $sequence = "";
	
	open(ORIGINAL, $infile) or die "Cannot open $infile\n";
	open(NEW, '>', $outfile) or die "Cannot open $outfile\n";


 	while ( my $line = <ORIGINAL> )
	{
		if ( $line =~ /^[a-zA-Z]/ ) # if line is sequence
		{
			chomp $line;
			$sequence = $sequence . $line ;
		}
		
		else #if allele doesn't exist in unique-hash
		{
			if ( !exists($unique{$sequence}) && length($sequence) > 1 ) # 
			{ 
    				print NEW ">" . $allelenumber . "\n" . $sequence . "\n"; 
				$unique{$sequence} = 1;
				$allelenumber++; 
			}
			$sequence = ""; 
		}
	}
 	close ORIGINAL; close NEW;
 	
 	$file =~ s{\.[^.]+$}{}; # remove extension from locus name
 	print UNIQUENUC $file . "," . ($allelenumber - 1) . "\n"; 

}

print "All done!\n";

#---------------------------------------------------------------
# Subroutines
#---------------------------------------------------------------

sub Usage( ; $ )
{
	my $message = $_[0] || '';

	print << 'EOU';

keepseenloci.pl

Description:
	  
Usage:
keepseenloci.pl [ options ]

Copies only the alleles in FASTA files that appear in a locus/isolate table. 

-din: Directory with FASTA files where "locusname.FAS" is the file name, like BACT000001.FAS
-dout: Directory where filtered FASTA files will be saved.
EOU

print "\n\nQuit because: $message\n\n";
print "ARGV was " . join (", ", @ARGV) . "\n";
}
