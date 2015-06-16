#!/usr/bin/perl

#--------------------------------------------------------
# PROGRAM: keepseenloci.pl
# AUTHOR:  Sofia Hauck
# CREATED: 09.06.2015
# UPDATED: ----------
# VERSION: v1.00
#--------------------------------------------------------
# VERSION HISTORY
# v1.00 (09.06.2015) created
#--------------------------------------------------------
# transposing adapted from http://stackoverflow.com/questions/3249508/transpose-in-perl

use strict;
use warnings;

# Declares subroutines
sub Usage( ; $ );

# Defines scalars needed from command line
my $fTab;
my $dFAS;

# Get Command line options, exits if conditions don't look right
if( scalar(@ARGV) < 1 ) { Usage("Not enough command line options"); exit; }
my $i = 0;
my $arg_cnt = 0; 
for ($i=0; $i<=$#ARGV; $i++)
{
	if($ARGV[$i] eq "-h")	       { Usage("You asked for help"); exit; }
	if($ARGV[$i] eq "-table")         { $fTab = $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-dfasta")        { $dFAS = $ARGV[$i+1] || ''; $arg_cnt++; }
}

# Command line option checks
if(! defined $fTab) { Usage("Missing Option: -table <FILE>"); exit; }
if(! defined $dFAS) { Usage("Missing Option: -dfasta DIRECTORY"); exit; }

# File checks
if(! -e $fTab)   { Usage("Input file does not exist: $fTab"); exit; }
if(! -e $dFAS)  { Usage("Output directory already exists: $dFAS"); exit; }


# Transposing:

my @original = (); #where data will go in the beginning
my @transposed = (); #where data is moved as its processed
my $columncount = 0; my $rowcount = 0;

# reads table into the @original array
open(INFILE, $fTab) or die "Cannot open $fTab\n";
	while ( my $line = <INFILE> )		
	{
        chomp($line); $line =~ s/\r\n?|\n//g; #just in case, removes all other types of line break
        
        my @row = split (',', $line);

        $rowcount = @row;
        $columncount ++;
        push ( @original, [ @row ] );
	} 
close(INFILE);

if ( $columncount == 1 )
{ print "You only have one column so your csv is probably not in Unix (LF) format.\n"; }

# actually does the transposition
for my $row (@original)
{
 	for my $column (0 .. $rowcount) 
 	{ push(@{$transposed[$column]}, $row->[$column]); }
}
@transposed = splice (@transposed, 0, $rowcount); #removes empty rows if more isolates than loci

# puts the transposed results somewhere sensible
my @results = (); 

for my $new_row (@transposed) 
{
	my $result; 
	for my $new_col (@{$new_row}) 
	{
		$result = $result . $new_col . ",";
	}
	$result = $result . "\n"; 
	push (@results, $result);
}


# Copying only relevant loci

open(INFILE, $fTab) or die "Cannot open $fTab\n";


foreach my $locusrow (@results)
{
	chomp($locusrow);
	my @row = split (',', $locusrow);

	# empty hash to find unique values, split line into elements of array, take out locus name
	my %unique_alleles = ();
	my @alleleArray = split (',', $locusrow);
	my $locusname = shift(@alleleArray);
	
	# makes a list of the unique alleles, including sorting them numerically
	foreach my $allele (@alleleArray) #go through each allele in locus
    	{
    		if ( !exists($unique_alleles{$allele})) #if allele doesn't exist in unique-hash
    		{ $unique_alleles{$allele} = 1; } #then add it with allele number as key, value as 1
	}

	my $originalFAS = $dFAS."/".$locusname.".FAS";
	my $reducedFAS = $dFAS."/reduced/".$locusname.".FAS";

	if ( open(FULLFASTA, $originalFAS) )
	{
		open(REDFASTA, '>', $reducedFAS) or die "Cannot open $reducedFAS\n";
		{       
			my $save = 1;
			while ( my $line = <FULLFASTA> )
			{
				if ( $line =~ /^>/ )
				{ 
					chomp $line; 
					my ($junk, $allelenumber) = split ('_', $line);
					if ( exists($unique_alleles{$allelenumber}) )
						{
							print REDFASTA $line, "\n";
							$save = 1;
						}
						elsif ( !exists($unique_alleles{$allelenumber}))
						{
							# print nothing
							$save = 0;
						}
				}
				elsif ( $line =~ /^[A-Z]/ )
				{
					if ($save == 1) { print REDFASTA $line; }
					elsif ($save == 0) { }
				}	
			}
		}  
		close(FULLFASTA);
		close(REDFASTA);
	}
	else { print "$locusname did not exist as a FASTA file.\n"; }	
	print "$locusname, ";
}
close(INFILE);


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

Results will be in folder with the FASTA files, in a subfolder called "reduced".
Actually you need to create this folder (and leave it empty) before running this script...

-table: Table as csv with rows as loci.
-dfasta: Directory with fasta files where "locusname.FAS" is the file name, like BACT000001.FAS

EOU

print "Quit because: $message\n";
print "ARGV was " . join (", ", @ARGV) . "\n";
}

