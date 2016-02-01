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



use strict;
use warnings;

# Declares subroutines
sub Usage( ; $ );
sub Unique( ); 

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

open(INGROUPS, $fGroup) or die "Cannot open $fGroup\n";
	while ( my $line = <INGROUPS> )
	{
		chomp($line); $line =~ s/\r\n?|\n//g; # clean up line breaks
		my @row = split ('\t', $line); # split row by tabs into array of cells
		$groups{$row[0]} = $row[1]; # uses only the first and second column, extra is ignored
	}
close(INGROUPS);

# check that only two groups exist exist
my @groupnames =  ( values %groups );
@groupnames = Unique ( @groupnames );


if ( @groupnames != 2 )
{ print "If there aren't just two groups, this isn't going to work.\n"; exit;}
else
{ print "The groups are '" . $groupnames[0] . "' and '" . $groupnames[1] . "'.\n"; }

# Read in the main table

# now that we know whether to transpose or not, read the table

my @aoaTable = ();
my $rowcount = 0;
	
	open(INTABLE, $fTable) or die "Cannot open $fTable\n";

	while ( my $line = <INTABLE> )
	{
		chomp($line); $line =~ s/\r\n?|\n//g; # clean up line breaks
		  
		my @row = split ('\t', $line); # split row by tabs into array of cells
		$rowcount ++;
		
		push ( @aoaTable, [ @row ] );
	} 

	if ( $rowcount <= 1 ) { die "You only have one column so something has probably gone wrong with line breaks.\n"; }


close(INTABLE);


# Determining which columns go to which group

my @headerrow = $aoaTable[0];

print @headerrow; 




# Actually the interesting bit of sorting alleles


# go through table that is now in AoA format
# each array 


# open(OUTFILE, '>>', $fOut); #open results file
# 	print OUTFILE "locus" ."\t". "status" ."\t". "total alleles" ."\t". "specific alleles" ."\t". "shared alleles" ."\t". "alleles shared"."\n";
# close(OUTFILE); 
# 
# # goes through each locus (as key in hash) and its array (list of allele numbers)
# my @hashkeys = keys(%locitable);
# @hashkeys = sort(@hashkeys); 
# 
# foreach my $locus ( @hashkeys )
# {
# 	# start with empty hash and array
# 	my %unique_alleles = ();
# 	my @sharedalleles = ();
# 	
# 	# take all the elements of the array (list of allelel numbers) at that key (locus)
# 	my @alleleArray = @{ $locitable{$locus} };
# 	
# 	# go through array for each element, making a frequency count hash 
# 	foreach my $allele (@alleleArray)
# 	{
#     		# create that allele number as a value in hash if it doesn't yet exist
#     		if ( !exists($unique_alleles{$allele}))
#     		{ $unique_alleles{$allele} = 1; }
#     		
#     		else
#     		{ 
#     			$unique_alleles{$allele} ++;
#     			if ($unique_alleles{$allele} == 2) 
#     			{ push @sharedalleles, $allele; } #when allele is known to be shared, add to @sharedallele
#     		} 
# 	}
# 	
# 	# take values (frequencies of alleles) from frequency count hash
# 	my @values = values(%unique_alleles);
# 	@values = sort {$a <=> $b} @values;
# 	
# 	# determine count of total, shared and specific alleles
# 	my $alleleCount = @values;
# 	my $specificCount = $alleleCount - @sharedalleles;
# 	my $sharedCount = @sharedalleles;
# 	my $status;
# 	
# 	# determining results
# 	
# 	# only one allele, means there is no variation and hence no specificity is possible
# 	if ($#values == 0)
# 	{ $status = "no variation"; }
# 	
# 	else
# 	{
# 		my $first = shift(@values); #first value in frequencies into $first
# 		my $last = pop(@values); #last value in frequencies into $last
# 		
# 		if ($first == $last) #if first and last are all the same, then all are the same since array is sorted
# 		{
# 			if ($first == 1) #if all alleles occur only once across all groups
# 			{ $status = "all specific"; } #then must be specific to the group
# 			elsif ($first == 2) #if all alleles occur in each group once
# 			{ $status = "all shared"; } #then must all overlap
# 		}
# 		else #mix of frequencies
# 		{
# 			if ($last == 2) #if highest frequency is total number of files compared
# 			{ $status = "shared"; }
# 			else #no alleles shared among all of the groups
# 			{ $status = "some shared"; }
# 		}
# 	}
# 	
# 	open(OUTFILE, '>>', $fOut); #open results file
# 	print OUTFILE $locus ."\t". $status ."\t". $alleleCount ."\t". $specificCount ."\t". $sharedCount ."\t". join(", ", @sharedalleles) ."\n";
# 	close(OUTFILE); 
# 	
# }


#---------------------------------------------------------------
# Subroutines
#---------------------------------------------------------------


sub Unique
 {
    my %seen;
    grep !$seen{$_}++, @_;
}

my @array = qw(one two three two three);
my @filtered = uniq(@array);

print "@filtered\n";

sub Usage( ; $ )
{
	my $message = $ARGV[0] || '';

	print << 'EOU';

allele-overlap.pl

Description:
	

Example Input Format:
  Dataset as exported from PubMLST or Genome Comparator
  
  Locus		Asdfg	Qwert	Zxcvb
  BACT000001	1		3		X
  BACT000002	I		new#1	0
  
  List of isolate names / ids and groups 
  Asdfg	L2
  Qwert	L2
  Zxcvb	L1
  
Example Output Format:
  
  
Usage:
Program.pl [ options ]

-intable        <FILE> - input filename
-ingroup        <FILE> - input filename
-out       <FILE> - output filename
-h|--help     - print usage instructions

EOU

print "Quit because: $message\n";
print "ARGV was..." . join ("\n", @ARGV);
}
