#!/usr/bin/perl

#--------------------------------------------------------
# PROGRAM: allele-overlap.pl
# AUTHOR:  Sofia Hauck
# CREATED: 20.06.2014
# UPDATED: 26.01.2016
# VERSION: v1.10
#--------------------------------------------------------
# VERSION HISTORY
# v1.00 (20.06.2014) created from spilicing older perl script together
# v1.10 (26.01.2016) making it work with big table and separating list 
#--------------------------------------------------------
# transposing adapted from http://stackoverflow.com/questions/3249508/transpose-in-perl
#--------------------------------------------------------

# start with a table and list that separates into two groups
# read the list into hash with keys as isolate names, values as name of list groups
	# check that only two values exist
# read the table, isolate by isolate
	# make an educated guess on which way around the table is?? ^[A-Z]{4}[0-9]{6}
	# either way, end up with array of arrays
	# for each array, take first item as locus name, others go as alleles into frequency hash
		# don't need to be integers
		# have specific checks for "X", "0" and "I"?
# only if isolate exists in grouping hash
# read into locus allele frequency hash for that group 
# do the usual results from there?
#
#







use strict;
use warnings;

# Declares subroutines
sub Usage( ; $ );

# Defines scalars needed from command line
my $fTable;
my $fGroup;
my $fOut;
my $skipSym;

# Get Command line options, exits if conditions don't look right
if( scalar(@ARGV) < 1 ) { Usage("Not enough command line options"); exit; }
my $i = 0;
my $arg_cnt = 0; 
for ($i=0; $i<=$#ARGV; $i++)
{
	if($ARGV[$i] eq "-h")			{ Usage("You asked for help"); exit; }
	if($ARGV[$i] eq "-table")		{ $fTable = $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-group")		{ $fGroup = $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-out")			{ $fOut = $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-skipsymbols")	{ $skipSym = 1 || ''; $arg_cnt++; }
}

# Check that required options were included
if(! defined $fTable) 	{ Usage("Missing Option: -table <FILE>"); exit; }
if(! defined $fGroup) 	{ Usage("Missing Option: -group <FILE>"); exit; }
if(! defined $fOut) 	{ Usage("Missing Option: -out <FILE>"); exit; }

# Check that files exist or don't exist as necessary
if(! -e $fTable)	{ Usage("Input file does not exist: $fTable"); exit; }
if(! -e $fGroup)	{ Usage("Input file does not exist: $fGroup"); exit; }
if(  -e $fOut )		{ Usage("Output file already exists: $fOut"); exit; }


# Transponsing if necessary, currently into array of arrays

if($transpose == 1) #!! change to detection?
{
	my @original = (); my @transposed = ();
	my $columncount = 0; my $rowcount = 0; # refers to rows and columns of the transposed table, not the original
	
	open(INTABLE, $fTable) or die "Cannot open $fTable\n";
		while ( my $line = <INTABLE> )
		{
		   chomp($line); $line =~ s/\r\n?|\n//g; # clean up line breaks
	   
		   my @row = split ('\t', $line); # split row by tabs into array of cells
		   if ( $rowcount < 2 ) { $rowcount = @row; } # transposed row count determined from header of original 
		
		   $columncount ++; #add to columncount with each row (will become column)
		   
		   push ( @original, [ @row ] );
		} 
	close(INTABLE);

	if ( $columncount == 1 ) { die "You only have one column so your csv is probably not in Unix (LF) format.\n"; }


	for my $row (@original) #actually does the transposition
	{
		for my $column (0 .. $rowcount) 
		{ push(@{$transposed[$column]}, $row->[$column]); }
	}
	
	@transposed = splice (@transposed, 0, $rowcount); #removes empty rows if more isolates than loci
}


#from this point on, if asks for infile, check if transpose, then use middle file?

open(INTABLE, $fTable) or die "Cannot open $fTable\n";

	<INTABLE>; # Skips header
	
	while ( my $line = <INTABLE> ) #goes through each row in the file at a time			
	{
		chomp($line); $line =~ s/\r\n?|\n//g; # clean up line breaks
        
		# somehow read line into array
		# my @linedata = grep {length()} split(/,/, $line); #put each allele into one item of an array, skips empty cells 

        	my $locus = shift(@linedata); #takes out first array item (locus name) puts it into $locus scalar
		    		
    		#for each allele, check if already exists in hash, creates if doesn't
  	  	my %unique_alleles = (); 
  	  	foreach my $allele (@linedata)
  	  	{
			if ( !exists($unique_alleles{$allele}))
				{ $unique_alleles{$allele} = 1; }
		}
				
		# Writes into outfile the name of the locus plus all its unique alleles 
		open(OUTFILE, '>>', $fOut) or die "Cannot open $fOut\n"; 
			print OUTFILE $locus;
			foreach my $key ( sort keys %unique_alleles ) 
			{
  				print OUTFILE "\t" . $key;
			}			
			print OUTFILE "\n";
		close(OUTFILE);		
	} 
close(INFILE);



#---------------------------------------------------------------
# Subroutines
#---------------------------------------------------------------

sub Usage( ; $ )
{
	my $message = $ARG[0] || '';

	print << 'EOU';

allele-overlap.pl

Description:
  Takes two lists, removes all items in second list from first list, returns remaining

Example Input Format:
  Dataset as exported from PubMLST
  
Example Output Format:
  		AB		BC		AC		ACB
  BACT01	unique	
  
Usage:
Program.pl [ options ]

--in        <FILE> - input filename
--subtract	<FILE>,<FILE>,... - any number of filenames, separated by commas
--out       <FILE> - output filename
-h|--help     - print usage instructions

EOU

print "Quit because: $message\n";
print "ARGV was " . join (", ", @ARGV);
}
