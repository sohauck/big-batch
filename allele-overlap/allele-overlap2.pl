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
my $transpose = "check needed"; 

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
	if($ARGV[$i] eq "-transpose")	{ $transpose = $ARGV[$i+1] || ''; $arg_cnt++; }

}

# Check that required options were included
if(! defined $fTable) 	{ Usage("Missing Option: -table <FILE>"); exit; }
if(! defined $fGroup) 	{ Usage("Missing Option: -group <FILE>"); exit; }
if(! defined $fOut) 	{ Usage("Missing Option: -out <FILE>"); exit; }

# Check that files exist or don't exist as necessary
if(! -e $fTable)	{ Usage("Input file does not exist: $fTable"); exit; }
if(! -e $fGroup)	{ Usage("Input file does not exist: $fGroup"); exit; }
if(  -e $fOut )		{ Usage("Output file already exists: $fOut"); exit; }


# Read in the "groups" file
# into a hash where isolates are keys, and groups are values
my %groups = ();

open(INGROUPS, $fTable) or die "Cannot open $fTable\n";
	while ( my $line = <INGROUPS> )
	{
	   chomp($line); $line =~ s/\r\n?|\n//g; # clean up line breaks

	   my @row = split ('\t', $line); # split row by tabs into array of cells

		$groups{@row[0]} = @row[1]; # uses only the first and second column, extra is ignored
	}
close(INGROUPS);

# check that only two values exist
my @groupnames = values %hash;
print "The groups are @groupsnames.\n";

if ( @groupnames != 2 )
{ print "If there aren't just two groups, this isn't going to work.\n"; exit;}


# Read in the main table

# if no manual override, read first line to check whether to transpose
if ( $transpose = "check needed" )
{
	print "Since direction of table isn't specific, now taking a guess...\n";
	
	open(INTABLE, $fTable) or die "Cannot open $fTable\n";
	
	my $header = <INTABLE>; 
	my @row = split ('\t', $line)
	
	if ( @row[0] = "Locus" ) # if the top right cell is "Locus" as Genome Comparator tables...
	{
		$transpose = 0;
		print "Since 'Locus' is at the top right, looks like GC format, and no transposing needed.\n";
	}
	elsif ( @row[0] = "id" ) # if the top right cell is "id" as most exported datasets are...
	{
		$transpose = 1;
		print "Since 'id' is at the top right, looks like 'Export Data set' format, and will be transposing.\n";
	}
	
	else # if the top right cell isn't anything helpful
	{
		my $itemcount = 0;
		my $lookslikelocus = 0;
		
		foreach @row # read each item in that first line
		{
			if ( /^[A-Z]{4}[0-9]{6}/ ) # see if it has the AAAA123456 format
			{ $lookslikelocus ++; }
			$itemcount ++;
		}
		
		if ( $lookslikelocus >= (0.9 * $itemcount) ) # if mostly matches to the locus name format
		{
			$transpose = 1;
			print "Looks like there are loci names in the header row, so will be transposing the table.\n";
		}
		else # if not very many matches for locus name format, the loci are probably in the rows already
		{
			$transpose = 0;
			print "Doesn't look like we have loci names in the header row, so no transposing needed.\n";
		}
	}
	
}

# now that we know whether to transpose or not, read the table

my @aoaTable = ();

if ( $transpose = 1 )
{
	my @original = (); my @aoaTable = ();
	my $columncount = 0; my $rowcount = 0; # refers to rows and columns of the transposed table, not the original
	
	while ( my $line = <INTABLE> )
	{
		chomp($line); $line =~ s/\r\n?|\n//g; # clean up line breaks
	   
		my @row = split ('\t', $line); # split row by tabs into array of cells
		if ( $rowcount < 2 ) { $rowcount = @row; } # transposed row count determined from header of original 
		
		$columncount ++; #add to columncount with each row (will become column)
		   
		push ( @original, [ @row ] );
	} 

	if ( $columncount == 1 ) { die "You only have one column so your csv is probably not in Unix (LF) format.\n"; }

	for my $row (@original) #actually does the transposition
	{
		for my $column (0 .. $rowcount) 
		{ push(@{$aoaTable[$column]}, $row->[$column]); }
	}
	
	@aoaTable = splice (@aoaTable, 0, $rowcount); #removes empty rows if more isolates than loci
}


if ( $transpose = 0 ) 
{
	my $rowcount = 0;
	
	while ( my $line = <INTABLE> )
	{
		chomp($line); $line =~ s/\r\n?|\n//g; # clean up line breaks
		  
		my @row = split ('\t', $line); # split row by tabs into array of cells
		$rowcount ++;
		
		push ( @aoaTable, [ @row ] );
	} 

	if ( $rowcount >= 1 ) { die "You only have no or just one column so your csv is probably not in Unix (LF) format.\n"; }
}

close(INTABLE);


# Determining which columns go to which group

# end up with two arrays with the indexes of the groups


# Actually the interesting bit of sorting alleles


# go through table that is now in AoA format
# each array 


open(OUTFILE, '>>', $fOut); #open results file
	print OUTFILE "locus" ."\t". "status" ."\t". "total alleles" ."\t". "specific alleles" ."\t". "shared alleles" ."\t". "alleles shared"."\n";
close(OUTFILE); 

# goes through each locus (as key in hash) and its array (list of allele numbers)
my @hashkeys = keys(%locitable);
@hashkeys = sort(@hashkeys); 

foreach my $locus ( @hashkeys )
{
	# start with empty hash and array
	my %unique_alleles = ();
	my @sharedalleles = ();
	
	# take all the elements of the array (list of allelel numbers) at that key (locus)
	my @alleleArray = @{ $locitable{$locus} };
	
	# go through array for each element, making a frequency count hash 
	foreach my $allele (@alleleArray)
	{
    		# create that allele number as a value in hash if it doesn't yet exist
    		if ( !exists($unique_alleles{$allele}))
    		{ $unique_alleles{$allele} = 1; }
    		
    		else
    		{ 
    			$unique_alleles{$allele} ++;
    			if ($unique_alleles{$allele} == 2) 
    			{ push @sharedalleles, $allele; } #when allele is known to be shared, add to @sharedallele
    		} 
	}
	
	# take values (frequencies of alleles) from frequency count hash
	my @values = values(%unique_alleles);
	@values = sort {$a <=> $b} @values;
	
	# determine count of total, shared and specific alleles
	my $alleleCount = @values;
	my $specificCount = $alleleCount - @sharedalleles;
	my $sharedCount = @sharedalleles;
	my $status;
	
	# determining results
	
	# only one allele, means there is no variation and hence no specificity is possible
	if ($#values == 0)
	{ $status = "no variation"; }
	
	else
	{
		my $first = shift(@values); #first value in frequencies into $first
		my $last = pop(@values); #last value in frequencies into $last
		
		if ($first == $last) #if first and last are all the same, then all are the same since array is sorted
		{
			if ($first == 1) #if all alleles occur only once across all groups
			{ $status = "all specific"; } #then must be specific to the group
			elsif ($first == 2) #if all alleles occur in each group once
			{ $status = "all shared"; } #then must all overlap
		}
		else #mix of frequencies
		{
			if ($last == 2) #if highest frequency is total number of files compared
			{ $status = "shared"; }
			else #no alleles shared among all of the groups
			{ $status = "some shared"; }
		}
	}
	
	open(OUTFILE, '>>', $fOut); #open results file
	print OUTFILE $locus ."\t". $status ."\t". $alleleCount ."\t". $specificCount ."\t". $sharedCount ."\t". join(", ", @sharedalleles) ."\n";
	close(OUTFILE); 
	
}


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
