#!/usr/bin/perl

#--------------------------------------------------------
# PROGRAM: align-varsites.pl
# AUTHOR:  Sofia Hauck
# CREATED: 09.06.2015
#--------------------------------------------------------

use strict;
use warnings;
$| = 1;

# Declares subroutines
sub Usage( ; $ );

# Get Command line options, exits if conditions don't look right
if( scalar(@ARGV) < 1 ) { Usage("Not enough command line options; need at least the input folder"); exit; }
my $dir = shift(@ARGV); # directory where all the FASTA files to be aligned are (can be nuc. or aa) 
my @mafftarg = @ARGV; # rest of @ARGV are arguments that will be passed directly to mafft

# Command line option & file checks
if(! defined $dir) { Usage("Missing Option: a directory that as the first argument"); exit; }
if(! -e $dir)   { Usage("Input directory does not exist: $dir"); exit; }


# get names of all the files in the directory
opendir (ORIGDIR, $dir) or die "Cannot open directory: $!";
	my @files = readdir ORIGDIR;
	@files = grep(/^([A-Z]|[a-z]|[0-9])/,@files);
	if ($#files < 1)
	{ die "It looks like you have no FASTA files to align. If their file names don't begin with letters or numbers they are being excluded.\n";}
closedir ORIGDIR;

# Loop that aligns all the files using MAFFT and puts them in new directory

# naming and making directory 
my @splitdir = split(/\//, $dir);
my $oldfolder = pop @splitdir;
my $alidir = join("/", @splitdir) . "/" . $oldfolder . "-a/";
if( -e $alidir)   { Usage("Output directory already exists: $alidir"); exit; }

mkdir $alidir;


# letting you know what's going to happen
print "Adding to MAFFT arguments: --clustalout --quiet for output in CLUSTAL format (necessary for non-variable count), and no STDOUT reports.\n";
push @mafftarg, ("--clustalout","--quiet");

print "Aligned up to...\n";

foreach my $file (@files)
{
 	my $infile = $dir . "/" . $file; # original file from folder given in command line
 	my $outfile = $alidir . $file; # aligned file into aligned folder, with same file name
 	
 	my $command = "mafft " . join (" ", @mafftarg) . " " . $infile . " > " . $outfile ;
 	system ($command); # passes mafft command to terminal
 	
 	# So you have something to watch while it runs... 
 	print "\r$file";
}
 

# Count of number of non-variable sites per locus
# essentially the number of asterisks per file

my @results = ();

foreach my $file (@files)
{
 	my $alignedfile = $alidir . $file;
	my $total = 0;
	
	# Loop per line inside loop per file, adds to asterisk count
	open (ALIGNED, $alignedfile) or die "Cannot open $alignedfile";
		    while (my $line = <ALIGNED>)
		    {
		    	my $count = ($line =~ tr/\*//);
		    	$total = $total + $count; 
		    }
	close ALIGNED;
	 
	my ($locusname, $extension) = split (/\./, $file);
	push @results, ($locusname . "," . $total);
}

# Then write results to a new file 
my $resultfile = substr($alidir, 0, -1) . "-count-nvs.txt";

open(RESULT, '>', $resultfile) or die "Cannot open $resultfile\n";
	print RESULT "locus,varsites\n";

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
Recommend "--maxiterate 0 --retree 1" for fastest option.

Output: 
Folder with same name and number of files but named -aligned with aligned FASTA files, and 
one output text file with the count of variable sites, in the parent directory of both previously mentioned folders.

EOU

print "Quit because: $message\n";
print "ARGV was " . join (", ", @ARGV) . "\n";
}

