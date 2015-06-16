#!/usr/bin/perl

#--------------------------------------------------------
# PROGRAM: align-varsites.pl
# AUTHOR:  Sofia Hauck
# CREATED: 09.06.2015
# UPDATED: ----------
# VERSION: v1.00
#--------------------------------------------------------
# VERSION HISTORY
# v1.00 (10.06.2015) created
#--------------------------------------------------------

use strict;
use warnings;

# Declares subroutines
sub Usage( ; $ );

# Defines options needed from command line


# Get Command line options, exits if conditions don't look right
if( scalar(@ARGV) < 2 ) { Usage("Not enough command line options; probably missing MAFFT arguments"); exit; }
my $dir = shift(@ARGV);
my @mafftarg = @ARGV; # rest of @ARGV are arguments that will be passed directly to mafft
print "MAFFT arguments are @mafftarg.\n";

# Command line option & file checks
if(! defined $dir) { Usage("Missing Option: a directory that as the first argument"); exit; }
if(! -e $dir)   { Usage("Input directory does not exist: $dir"); exit; }

# get names of all the files in the directory
opendir (ORIGDIR, $dir) or die "Cannot open directory: $!";
	my @files = readdir ORIGDIR;
	@files = grep(/^([A-Z]|[a-z]|[0-9])/,@files);
	if ($#files < 1)
	{ die "It looks like you have no FASTA files to align. If their file names don't begin with letters or numbers they are being removed.\n";}
closedir ORIGDIR;

# Loop that aligns all the files using MAFFT and puts them in new directory

# naming and making directory 
my @splitdir = split(/\//, $dir);
my $oldfolder = pop @splitdir;
my $alidir = join("/", @splitdir) . "/" . substr($oldfolder, 0, 5) . "-aligned/";
mkdir $alidir;

# letting you know what's going to happen
print "Adding to MAFFT arguments: output in CLUSTAL format (for non-variable count), and quiet terminal output.\n";
push @mafftarg, ("--clustalout", "--quiet");

print "Aligned up to...";

foreach my $file (@files)
{
 	my $infile = $dir . "/" . $file; # original file from folder given in command line
 	my $outfile = $alidir . $file; # aligned file into aligned folder, with same file name
 	
 	my $command = "mafft " . join (" ", @mafftarg) . " " . $infile . " > " . $outfile ;
 	system ($command); # passes mafft command to terminal
 	
 	# So you have something to watch while it runs... 
 	print $file . ", ";
}
 

# Count of number of non-variable nucleotides per locus
# essentially the number of askerisks per file

my @results = ();

foreach my $file (@files)
{
 	my $alignedfile = $alidir . $file;
	my $totalcount = 0;
	
	# Loop per line inside loop per file, adds to asterisk count
	open (ALIGNED, $alignedfile) or die "Cannot open $alignedfile";
		    while (my $line = <ALIGNED>)
		    {
		    	my $count = ($line =~ tr/\*//);
		    	$totalcount = $totalcount + $count; 
		    }
	close ALIGNED;
	
	push @results, ($file . "\t" . $totalcount);
}

# Then write results to a new file 
my $resultfile = substr($alidir, 0, -1) . "var-count.txt";

open(RESULT, '>', $resultfile) or die "Cannot open $resultfile\n";
	print RESULT "Count of number of variable sites per locus in $alidir files.\n\n";

	foreach my $result (@results)
	{ print RESULT $result , "\n"; } 
close RESULT;

#---------------------------------------------------------------
# Subroutines
#---------------------------------------------------------------

sub Usage( ; $ )
{
	my $message = $_[0] || '';

	print << 'EOU';

align-varsites.pl
follows from keenseenloci.pl

Description:
	  
Usage:
align-varsites.pl [ options ]

Takes FASTA files from a folder, aligns them with MAFFT then saves results.
Then goes back over files and checks for number of non-variable loci per aligment.

Input: 
First argument is folder with FASTA files, all others are passed on to MAFFT.

Output: 
Folder with same name and number of files but named -aligned with aligned FASTA files, and 
one output text file with the count of variable sites, in the parent directory of both previously mentioned folders.

EOU

print "Quit because: $message\n";
print "ARGV was " . join (", ", @ARGV) . "\n";
}

