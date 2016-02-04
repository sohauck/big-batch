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
	

use strict;
use warnings;

# Declares subroutines
sub Usage( ; $ );
# also Unique, ReadTableIn, Transpose, etc...

# Defines scalars needed from command line
my $fTable;
my $fGroup;
my $fOut;
my $skipSym = 1;
my @symbols = ( "X", "I", "0" ); 
my $transpose = "check needed"; # unless manually specified, will take a guess

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
		
		# first column is isolate unique ID, second column is whatever group, ignore the rest
		$groups{$row[0]} = $row[1]; 
	}
close(INGROUPS);

# check that only two groups exist exist
my @groupnames =  ( values %groups ); # groups are the values in the hash that had isolate-group pairs
@groupnames = Unique ( @groupnames ); # only really care about the count of 2


if ( @groupnames != 2 )
{ print "If there aren't just two groups, this isn't going to work.\n"; exit;}
else
{ print "The groups are '" . $groupnames[0] . "' and '" . $groupnames[1] . "'.\n"; }




# Checking if transposing is needed if nothing is specified 

if ( $transpose eq "check needed" )
{
	print "Since direction of table was not specified, now taking a guess...\n";
	
	# take just the header row from the table file
	open(INTABLE, $fTable) or die "Cannot open $fTable\n";
	my $header = <INTABLE>; close(INTABLE);
	my @headerrow = split ('\t', $header);
	
	if ( $headerrow[0] eq "Locus" ) # if the top right cell is "Locus" as Genome Comparator tables...
	{
		$transpose = 0;
		print "Since 'Locus' is at the top right, looks like GC format, and no transposing needed.\n";
	}
	elsif ( $headerrow[0] eq "id" ) # if the top right cell is "id" as most exported datasets are...
	{
		$transpose = 1;
		print "Since 'id' is at the top right, looks like 'Export Data set' format, and will be transposing.\n";
	}
	
	else # if the top right cell isn't anything helpful
	{
		my $itemcount = 0;
		my $lookslikelocus = 0;
		
		foreach ( @headerrow ) # read each item in that first line
		{
			if ( /^[A-Z]{4}[0-9]{6}/ ) # see if it has the AAAA123456 format
			{ $lookslikelocus ++; }
			$itemcount ++;
		}
		
		if ( $lookslikelocus >= (0.75 * $itemcount) ) # if mostly matches to the locus name format (to account for fragments)
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



# Read in the main table

my $tableref = ReadTableIn ( $fTable );
my @aoaTable = @$tableref; #puts table, in array of arrays format, into aoaTable
if ( $transpose == 1 )
{
	$tableref = TransposeTable ( \@aoaTable );
	@aoaTable = @$tableref;
}

# Determining which columns go to which group

my $headerref = $aoaTable[0]; #dereferences the first line of the table
my @grouporder = (); # @grouporder should contain only the two names that exist in @groupnames, and undef otherwise

foreach ( @$headerref ) # loops through each item in first row
{
	if ( /^(\S+)/ && $groups{$1} ) # takes only the first "word", and only if it exists in the %groups hash
	{ push @grouporder, $groups{$1}; } # silently add name of group to array as in order
	else 
	{ 
		push @grouporder, "Neither Group";
		print "'" . $_ . "' was not an isolate found in the grouping list.\n"; 
	}
}


# setting up the results file
open(OUTFILE, '>>', $fOut); #open results file
	print OUTFILE 	"Locus" 				."\t". 
					"Status" 				."\t". 
					"Total allele count" 	."\t".
					"Shared allele count" 	."\t". 
					"Shared alleles" 		."\t". 
					"Alleles specific to $groupnames[0]"	."\t". 
					"Alleles specific to $groupnames[1]"	."\n";


# Actually the interesting bit of sorting alleles

for (my $i = 1; $i <= $#aoaTable; $i++) # going through each locus (one per line) in input table, skips header row (i=0)
{
	my $lineref = $aoaTable[$i]; # find reference to line (array) in AoA table
		
	my $locus = @$lineref[0]; # put first item in line into $locus scalar, since it is the locus ID
	$locus =~ /^(\S+)/; $locus = $1; # keep only what's before the first white space (in case aliases etc. are there too)
	
	my %alleleObs = (); # start with an empty hash, will have allele ID as keys, 1 / 2 / 3 as values
	
	for (my $i = 0; $i <= $#grouporder; $i++) # going through each item in line (usually alleles)
	{
		my @alleles = split (/,|;/, @$lineref[$i]); # in case of paralogous loci
	
		if ( $skipSym == 1 )
		{
			my %symbols; @symbols{@symbols} = ();
			@alleles = grep !exists $symbols{$_}, @alleles;
		}
		
		if ( scalar(@alleles) == 0 ) { next; } # don't bother with rest of loop if there are no alleles
				
		if ( $grouporder[$i] eq $groupnames[0] ) # if column is from isolate in group 1
		{ 			
			foreach ( @alleles )
			{
				if ( !exists $alleleObs{$_} ) # if doesn't exist at all, add as existing in group 1
				{ $alleleObs{$_} = 1; }
				elsif ( $alleleObs{$_} == 2 ) # if exists only as group 2 so far, add as existing in both
				{ $alleleObs{$_} = 3; }
			}
		}
		
		elsif ( $grouporder[$i] eq $groupnames[1] ) # if column is from isolate in group 2
		{ 
			foreach ( @alleles )
			{
				if ( !exists $alleleObs{$_} ) # if doesn't exist at all, add as existing in group 2
				{ $alleleObs{$_} = 2; }
				elsif ( $alleleObs{$_} == 1 ) # if exists only as group 1 so far, add as existing in both
				{ $alleleObs{$_} = 3; }
			}
		}
	} # closes item in line loop
	
	my @only1 = (); my @only2 = (); my @both = (); 

	foreach ( keys %alleleObs )
 	{		
		if ( $alleleObs{$_} == 1 )
		{ push @only1, $_ }
 		if ( $alleleObs{$_} == 2 )
 		{ push @only2, $_ }
 		if ( $alleleObs{$_} == 3 )
 		{ push @both , $_ }
	}
	
	# determining status
	my $status; 
	
	if ( (@only1+@only2+@both) <= 1 )
	{ $status = "No variation"; }
	elsif ( (@only1+@only2) == 0 )
	{ $status = "All shared"; }
	elsif ( @both == 0 )
	{ $status = "Specific"; }
	else 
	{ $status = "Mixed"; }
		
	print OUTFILE 	$locus					."\t". # Locus
					$status 				."\t". # Status
					(@only1+@only2+@both) 	."\t". # Total allele count
					scalar(@both)	 		."\t". # Shared allele count
					join (";" , @both )		."\t". # Shared alleles
					join (";" , @only1 )	."\t". # Alleles specific to $groupnames[0]
					join (";" , @only2 )	."\n"; # Alleles specific to $groupnames[1]

}  # closes locus loop


close(OUTFILE); 



#---------------------------------------------------------------
# Subroutines
#---------------------------------------------------------------


sub Unique
 {
    my %seen;
    grep !$seen{$_}++, @_;
}

sub ReadTableIn
{
	my $infile = $_[0];
	my @aoaTable = ();
	my $rowcount = 0; 
	
	open(INTABLE, $infile) or die "Cannot open $infile\n";

		while ( my $line = <INTABLE> )
		{
			chomp($line); $line =~ s/\r\n?|\n//g; # clean up line breaks
		  
			my @row = split ('\t', $line); # split row by tabs into array of cells
			$rowcount ++;
		
			push ( @aoaTable, [ @row ] );
		} 

		if ( $rowcount <= 1 ) { die "You only have one column so something has probably gone wrong with line breaks.\n"; }
	close(INTABLE);
	
	return \@aoaTable;
}

sub TransposeTable
{	
	print "Passed $_[0]   ";
	
	my $aoaRef = $_[0]; 
	my @original = @$aoaRef; # @original is now array of references to arrays
	
	print "Original now holds @original";
		
	my @transposed = ();
	my $columncount = 0; 
		
	for my $row (@original) 
	{
  		for my $column (0 .. $#{$row}) 
  		{
    		push(@{$transposed[$column]}, $row->[$column]);
    		$columncount = $#{$row};
  		}
	}
	
	@transposed = splice (@transposed, 0, $columncount); #removes empty rows if more isolates than loci
	
	return \@transposed; 
}

sub Usage( ; $ )
{
	my $message = $ARGV[0] || '';

	print << 'EOU';

allele-overlap.pl

Description:
	
	Symbols that are excluded: X (missing), I (incomplete), 0 (zero)

Example Input Format:
  Dataset as exported from PubMLST or Genome Comparator
  Table of isolates x loci 
  
  Locus		Asdfg	Qwert	Zxcvb
  BACT000001	1		3		X
  BACT000002	I		new#1	0
  
  List of isolate names / ids and groups 
  columns beyond the first and second will be ignored
  
  Asdfg	L2
  Qwert	L2
  Zxcvb	L1
  
  
Example Output Format:
  (not too sure yet...)
  
Usage:
Program.pl [ options ]

-intable        <FILE> - input filename
-ingroup        <FILE> - input filename
-out       <FILE> - output filename
-h|--help     - print usage instructions

EOU

print "Quit because: $message\n";
print "ARGV was..." . join ("; ", @ARGV);
}
