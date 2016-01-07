#!/usr/bin/perl

#--------------------------------------------------------
# PROGRAM: keepseenloci.pl
# AUTHOR:  Sofia Hauck
# CREATED: 09.06.2015
#--------------------------------------------------------
# transposing adapted from http://stackoverflow.com/questions/3249508/transpose-in-perl

use strict;
use warnings;
$| = 1;

# Declares subroutines
sub Usage( ; $ );

# Defines scalars needed from command line
my $fTab;    # file with isolate/locus table where relevant loci are named
my $dFAS;    # directory where all the FASTA files for each locus are
my $dOut;    # directory where only the relevant parts of those FASTA files will be copied to
my $dup = 1; # turned to the smallest number of duplicates necessary for locus to be considered, default is 1

# Get Command line options, exits if conditions don't look right
if( scalar(@ARGV) < 1 ) { Usage("Not enough command line options"); exit; }
my $i = 0;
my $arg_cnt = 0; 
for ($i=0; $i<=$#ARGV; $i++)
{
	if($ARGV[$i] eq "-h")	    { Usage("You asked for help"); exit; }
	if($ARGV[$i] eq "-tin")           { $fTab = $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-din")           { $dFAS = $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-dout")          { $dOut = $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-dup")           { $dup = $ARGV[$i+1] || ''; $arg_cnt++; }

}

# Command line option checks & File checks
if(! defined $fTab)  { Usage("Missing Option: -table <FILE>"); exit; }
if(! defined $dFAS)  { Usage("Missing Option: -dfasta DIRECTORY"); exit; }
if(! defined $dOut)  { Usage("Missing Option: -out DIRECTORY"); exit; }

if(! -e $fTab)  { Usage("Input file does not exist: $fTab"); exit; }
if(! -e $dFAS)  { Usage("Input directory doesn't exist: $dFAS"); exit; }
if(  -e $dOut)  { Usage("Output directory already exists: $dOut"); exit; }

# Makes output directory and prepares output count file, including header
mkdir $dOut;
my $uniquefile = $dOut . "-count-nuc.txt";

if( -e $uniquefile)  { Usage("Output file already exists: $uniquefile"); exit; }

open(UNIQUENUC, '>', $uniquefile) or die "Cannot open $uniquefile\n";
	print UNIQUENUC "locus,count-nuc,missing\n";


# Transposing, starting with set up...
print "\nTransposing...";
my @original = (); #where data will go in the beginning
my @transposed = (); #where data is moved as its processed
my $columncount = 0; my $rowcount = 0;

# reads table into the @original array
open(INFILE, $fTab) or die "Cannot open $fTab\n";
	while ( my $line = <INFILE> )		
	{
		chomp($line); $line =~ s/\r\n?|\n//g; #just in case, removes all other types of line break
	
		my @row = split (/[,\t]/, $line); # works if comma-separated file

		$rowcount = @row;
		$columncount ++;
		push ( @original, [ @row ] );
	} 
close(INFILE);

if ( $columncount == 1 ) { die "You only have one column so your csv is probably not in Unix (LF) format.\n"; }

# actually does the transposition
for my $row (@original)
{
 	for my $column (0 .. $rowcount) 
 	{ push(@{$transposed[$column]}, $row->[$column]); }
}
@transposed = splice (@transposed, 0, $rowcount); #removes empty rows if more isolates than loci

# puts the transposed results into an array of arrays
my @newtable = ();

for my $new_row (@transposed)
{
	my @result;
	for my $new_col (@{$new_row}) 
	{ 
		if ( !defined($new_col) ) # if empty cells for missing allele
		{ push (@result, "0"); } # write as "0", creates less confusion with undefined or empty variables
		
		else
		{ push (@result, $new_col); }
	}
	push (@newtable, \@result); # puts reference to this array into array of arrays
}

print " complete!\n"; # Transposing is done


# Copying only relevant loci
print "Now filter-copying...\n";

foreach my $locusrow (@newtable) # loop per locus
{
	# empty hash to find unique values, set locus name
	my %unique_alleles = ();
	
	my $locusname = shift(@{$locusrow});
		
	# populates that hash with unique allele numbers as keys and their number of appearances as values
	foreach my $allele (@{$locusrow}) #go through each allele in locus
    	{
    		if ( $allele =~ /;/ ) # in case of paralogous loci
    		{  
    			my @paralogous = split (/;/, $allele);
    			    			
    			foreach my $paraallele (@paralogous) #go through each allele in locus
			{ $unique_alleles{$paraallele}++; } #increase the frequency count	
    		}
    		
    		elsif ( $allele !~ /^\d+$/ ) # in case it isn't a numeric allele designation
    		{ $unique_alleles{"0"}++; } # just add as 0
    		
    		else # if just one numeric value, add it 
    		{ $unique_alleles{$allele}++; } 
	}
	
	my $originalFAS = $dFAS."/".$locusname.".FAS";
	my $reducedFAS = $dOut."/".$locusname.".FAS";

	if ( open(FULLFASTA, $originalFAS) )
	{
		open(REDFASTA, '>', $reducedFAS) or die "Cannot open $reducedFAS\n";
		{       
			my $save = 1; # whether to copy the sequence following a >identifier line
			my $count = 0; # counts how many alleles are copied over, for count-nuc file
			while ( my $line = <FULLFASTA> )
			{
				if ( $line =~ /^>/ ) # if line is >identifier
				{
					chomp $line; 
					my @alleletitle = split ('_', $line); #if in BIGS BACT01_1 format
					my $allelenumber = pop @alleletitle;
					$allelenumber =~ s/[>]+//g; #in case > is still in there if there is no other ID
					
					if ( exists($unique_alleles{$allelenumber}) ) # exists in hash of wanted loci
					{

						if ($unique_alleles{$allelenumber} >= $dup ) # and in frequence at or above cutoff
						{
							print REDFASTA "\n", $line, "\n";
							$save = 1;
							$count ++;
							$unique_alleles{$allelenumber} = 0; # set frequency to 0 as check that was copied
						}
						else 
						{ $save = 0; } # knows to skip the sequences lines that follow unwanted identifiers
					}
					else { $save = 0; }
				}
				
				elsif ( $line =~ /^[A-Z]/ ) # if is a sequence line, copy only is "save" is turned on by wanted identifier
				{
					chomp $line; 
					if ($save == 1)    { print REDFASTA $line; }
					elsif ($save == 0) { }
				}	
			}
			
			my $missing = 0; # sets the missing count back to empty
			
			if ( exists($unique_alleles{"0"}) ) # move count away from empty if any missing alleles were seen
			{ $missing = $unique_alleles{"0"}; } # gives the value in the frequency hash when the key is allele "0", the missing allele
			print UNIQUENUC $locusname . "," . $count . "," . $missing . "\n"; 

		}  
		close(FULLFASTA);
		close(REDFASTA);
	
		# check if all alleles seen in table were copied 
		foreach my $key ( sort keys %unique_alleles ) 
		{ 
			if ( $unique_alleles{$key} >= $dup && $key != 0 ) # if frequency is still above cut-off (copied ones are reset to 0) and isn't the missing allele, 0
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

-tin: Table as csv with columns as loci.
-din: Directory with FASTA files where "locusname.FAS" is the file name, like BACT000001.FAS
-dout: Directory where filtered FASTA files will be saved.
-dup: Can set a mininum frequence for allele appearance (default is 1) 
EOU

print "\n\nQuit because: $message\n\n";
print "ARGV was " . join (", ", @ARGV) . "\n";
}
