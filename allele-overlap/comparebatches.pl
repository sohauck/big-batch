#!/usr/bin/perl

#--------------------------------------------------------
# PROGRAM: comparebatches.pl
# AUTHOR:  Sofia Hauck
# CREATED: 16.09.2015
#--------------------------------------------------------

use strict;
use warnings;
$| = 1;

# Declares subroutines
sub Usage( ; $ );

# Defines scalars needed from command line
my $dIn1;    # directory where the FASTA files for each locus are, the one whose loci will be used
my $dIn2;    # directory where the FASTA files for each locus are,
my $dOut;    # directory where only the relevant parts of those FASTA files will be copied to

# Get Command line options, exits if conditions don't look right
if( scalar(@ARGV) < 1 ) { Usage("Not enough command line options"); exit; }
my $i = 0;
my $arg_cnt = 0; 
for ($i=0; $i<=$#ARGV; $i++)
{
	if($ARGV[$i] eq "-h")	          { Usage("You asked for help"); exit; }
	if($ARGV[$i] eq "-din1")          { $dIn1 = $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-din2")          { $dIn2 = $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-dout")          { $dOut = $ARGV[$i+1] || ''; $arg_cnt++; }
}

# Command line option checks & File checks
if(! defined $dIn1)  { Usage("Missing Option: -din1 DIRECTORY"); exit; }
if(! defined $dIn2)  { Usage("Missing Option: -din2 DIRECTORY"); exit; }
if(! defined $dOut)  { Usage("Missing Option: -dout DIRECTORY"); exit; }

if(! -e $dIn1)  { Usage("Input directory doesn't exist: $dIn1"); exit; }
if(! -e $dIn2)  { Usage("Input directory doesn't exist: $dIn2"); exit; }
if(  -e $dOut)  { Usage("Output directory already exists: $dOut"); exit; }

# Makes output directory and prepare output summary file, including header
mkdir $dOut;

my $uniquefile = $dOut . "-compare.txt";

if( -e $uniquefile)  { Usage("Output file already exists: $uniquefile"); exit; }

open(RESULTS, '>', $uniquefile) or die "Cannot open $uniquefile\n";
	print RESULTS "locus,count-1,count-2,both,1-only,2-only\n";

my @results = ();

# Get names of all the files in the directory
opendir (FIRSTFOLDER, $dIn1) or die "Cannot open directory: $!";
	my @files = readdir FIRSTFOLDER;
	@files = grep(/^([A-Z]|[a-z]|[0-9])/,@files); # in case of hidden files or anything funky
	if ($#files < 1)
	{ die "It looks like you have no FASTA files to align. If their file names don't begin with letters or numbers they are not being included.\n";}
closedir FIRSTFOLDER;


# Loop for each file that exists in the first directory, compares existence in both folders

foreach my $file (@files)
{
 	print "\r$file"; # so you have something to watch while it runs
 	
 	my $infile1 = $dIn1  . "/" . $file; # original file from folder given in command line
 	my $infile2 = $dIn2  . "/" . $file; # original file from folder given in command line
 	my $outfile = $dOut . "/" . $file; # new file into output folder, with only exceptions
 	
	my %unique = (); # empty hash to search for uniqueness 
	my $allelecount = 0; # empty scalar to count how many alleles are read in
	
	$file =~ s{\.[^.]+$}{};
	my $resultline = $file . ","; # locus name at the start of this line in results file
	
	open(FIRST, $infile1) or die "Cannot open $infile1\n";
 	while ( my $line = <FIRST> )
	{
		if ( $line =~ /^[a-zA-Z]/ ) # if line is sequence
		{
			$unique{$line} = 1;
			$allelecount++;
		}
	}
	close FIRST;
	$resultline = $resultline . $allelecount . ","; # add value for count-1 column
	
	$allelecount = 0;
	
	if ( open(SECOND, $infile2) )
	{
	 	while ( my $line = <SECOND> )
		{
			if ( $line =~ /^[a-zA-Z]/ ) # if line is sequence
			{
				if ( !exists($unique{$line}) && length($line) > 1 )
				{ $unique{$line} = 2; }
				
				elsif ( exists($unique{$line}) ) 
				{ $unique{$line} = 3; }
			
				$allelecount++;
			}
		}
		close SECOND;
		$resultline = $resultline . $allelecount . ","; # add value for count-2 column

	}
	
	else 
	{
		print "Cannot open $infile2, skipped\n";
		$resultline = $resultline . "0,"; # add value for count-2 column
	}
	
		
	# go through 
	my $countboth  = 0; # ==3
	my $count1only = 0; # ==1
	my $count2only = 0; # ==2
	
	open(UNMATCHED, '>', $outfile) or die "Cannot open $outfile\n";
	
	print UNMATCHED "Hello there!\n";
	
	foreach my $key ( sort keys %unique ) 
	{ 
		print "The value is " . $unique{$key};
		
		if ( $unique{$key} == 3 )
		{ $countboth ++; }
		if ( $unique{$key} == 1 )
		{ 	
			$count1only ++;
			print UNMATCHED ">FIRST_" . $count1only . "\n" . $key . "\n"; 
		}
		if ( $unique{$key} == 2 )
		{ 	$count2only ++;
			print UNMATCHED ">SECOND_" . $count2only . "\n" . $key . "\n";
		}
		
		
	}			
	
	$resultline = $resultline . $countboth . ",". $count1only . "," . $count2only; # add values for columns
	
	push ( @results, ($resultline)); 
	
	close UNMATCHED;

}

# sort results by locus name (as first item in row)
# print all results to results file
	@results = sort {$a cmp $b} @results;
	print RESULTS join("\n", @results);
close RESULTS; 


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
comparebatches.pl [ options ]

Compares FASTA files for loci in two folders, to check for overlap and differences. 

-din1: Directory with FASTA files where "locusname.FAS" is the file name, like BACT000001.FAS
-din2: likewise for din1 but the second folder.
-dout: Directory where exception FASTA files will be saved. Basis for results text file. 
EOU

print "\n\nQuit because: $message\n\n";
print "ARGV was " . join (", ", @ARGV) . "\n";
}


