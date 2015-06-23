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
$| = 1;

# Declares subroutines
sub Usage( ; $ );

# Defines scalars needed from command line
my $fTab; # file with isolate/locus table where relevant loci are named
my $dFAS; # directory where all the FASTA files for each locus are
my $dOut; # directory where only the relevant parts of those FASTA files will be copied to
my $dup = 1; #  turned to the smallest number of duplicates necessary for locus to be considered, default is 1

# Get Command line options, exits if conditions don't look right
if( scalar(@ARGV) < 1 ) { Usage("Not enough command line options"); exit; }
my $i = 0;
my $arg_cnt = 0; 
for ($i=0; $i<=$#ARGV; $i++)
{
	if($ARGV[$i] eq "-h")	          { Usage("You asked for help"); exit; }
	if($ARGV[$i] eq "-in")            { $fTab = $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-dfasta")        { $dFAS = $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-out")           { $dOut = $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-dup")           { $dup = $ARGV[$i+1] || ''; $arg_cnt++; }

}

# Command line option checks
if(! defined $fTab) { Usage("Missing Option: -table <FILE>"); exit; }
if(! defined $dFAS) { Usage("Missing Option: -dfasta DIRECTORY"); exit; }
if(! defined $dOut) { Usage("Missing Option: -out DIRECTORY"); exit; }

# File checks
if(! -e $fTab)		{ Usage("Input file does not exist: $fTab"); exit; }
if(! -e $dFAS)		{ Usage("Input directory doesn't exist: $dFAS"); exit; }
if(  -e $dOut)		{ Usage("Output directory already exists: $dOut"); exit; }

mkdir $dOut;

# Transposing:
print "\nTransposing...";
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
my @newtable = (); 

for my $new_row (@transposed) 
{
	my @result; 
	for my $new_col (@{$new_row}) 
	{
		push (@result, $new_col);
	}
	my $newline = join (",", @result);
	push (@newtable, $newline);
}

print " complete!\n";


# Copying only relevant loci

# naming and making count file
my $uniquefile = $dOut . "-count-nuc.txt";
if( -e $uniquefile)   { Usage("Output file already exists: $uniquefile"); exit; }
open(UNIQUENUC, '>', $uniquefile) or die "Cannot open $uniquefile\n";
	print UNIQUENUC "locus,count-nuc\n";


open(INFILE, $fTab) or die "Cannot open $fTab\n";

print "Now filter-copying...";
foreach my $locusrow (@newtable)
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
    		# in case of paralogous loci
    		if ( $allele =~ /;/ )
    		{
    			my @paralogous = split (';', $allele);
    			
    			foreach my $paraallele (@paralogous) #go through each allele in locus
			{
				if ( !exists($unique_alleles{$paraallele})) #if allele doesn't exist in unique-hash
    				{ $unique_alleles{$paraallele} = 1; } #then add it with allele number as key, value as 1
			}
    		}
    		
    		if ( !exists($unique_alleles{$allele})) #if allele doesn't exist in unique-hash
    		{ $unique_alleles{$allele} = 1; } #then add it with allele number as key, value as 1
    		
    		if (  exists($unique_alleles{$allele})) #if allele is already in the hash
    		{ $unique_alleles{$allele}++; } #increase the frequency count

	}

	my $originalFAS = $dFAS."/".$locusname.".FAS";
	my $reducedFAS = $dOut."/".$locusname.".FAS";

	if ( open(FULLFASTA, $originalFAS) )
	{
		open(REDFASTA, '>', $reducedFAS) or die "Cannot open $reducedFAS\n";
		{       
			my $save = 1;
			my $count;
			while ( my $line = <FULLFASTA> )
			{
				if ( $line =~ /^>/ )
				{
					chomp $line; 
					my ($locusname, $allelenumber) = split ('_', $line);
					if ( exists($unique_alleles{$allelenumber}) )
						{
							print REDFASTA "\n", $line, "\n";
							$unique_alleles{$allelenumber} = 0;
							$save = 1;
							$count ++;
						}
						elsif ( !exists($unique_alleles{$allelenumber}))
						{
							# print nothing
							$save = 0;
						}
				}
				elsif ( $line =~ /^[A-Z]/ )
				{
					chomp $line; 
					if ($save == 1)    { print REDFASTA $line; }
					elsif ($save == 0) { }
				}	
			}
			 
			print UNIQUENUC $locusname . "," . $count . "\n";

		}  
		close(FULLFASTA);
		close(REDFASTA);
	
		# check if all alleles were used
		foreach my $key ( sort keys %unique_alleles ) 
		{ 
			if ( $unique_alleles{$key} == 1 && $key != 0 )
			{ print "Did not find sequence for locus $locusname, allele $key.\n"; }
		}			
	}
	
	else { # if there isn't a FASTA file to copy over, give warning plus example alleles (might be just "0" which would explain missing file)
		print "$locusname did not exist as a FASTA file. Alleles in table were..."; 
		my $count2 = 0; 
		foreach my $key ( sort keys %unique_alleles ) 
		{ 
			if ( $count2 < 5 )
			{ print " $key"; $count2++; }
			else # if already printed 5 allele numbers, just leave it
			{ print " etc."; last; }
		 }
		print "\n";
	}	

} # closes per-locus loop


close(INFILE);
close(UNIQUENUC);
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

-in: Table as csv with rows as loci.
-dfasta: Directory with FASTA files where "locusname.FAS" is the file name, like BACT000001.FAS
-out: Directory where filtered FASTA files will be saved.

EOU

print "Quit because: $message\n";
print "ARGV was " . join (", ", @ARGV) . "\n";
}

