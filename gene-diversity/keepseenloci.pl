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
sub FreqCount();

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
	if($ARGV[$i] eq "-intab")            { $fTab = $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-din")        { $dFAS = $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-dout")           { $dOut = $ARGV[$i+1] || ''; $arg_cnt++; }
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

print "Now filter-copying...\n";
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
    		# if just one value
    		$unique_alleles{$allele}++; #increase the frequency count
    		
    		# in case of paralogous loci
    		if ( $allele =~ /;/ )
    		{
    			my @paralogous = split (';', $allele);
    			
    			foreach my $paraallele (@paralogous) #go through each allele in locus
			{ $unique_alleles{$allele}++; } #increase the frequency count
			
    		}

	}

	my $originalFAS = $dFAS."/".$locusname.".FAS";
	my $reducedFAS = $dOut."/".$locusname.".FAS";

	if ( open(FULLFASTA, $originalFAS) )
	{
		open(REDFASTA, '>', $reducedFAS) or die "Cannot open $reducedFAS\n";
		{       
			my $save = 1; # whether to copy the sequence following a >identifier line
			my $count; # counts how many alleles are copied over, for count-nuc file
			while ( my $line = <FULLFASTA> )
			{
				if ( $line =~ /^>/ ) # if line is >identifier
				{
					chomp $line; 
					my ($locusname, $allelenumber) = split ('_', $line);
					
					if ( exists($unique_alleles{$allelenumber}) )
					{
						if ($unique_alleles{$allelenumber} >= $dup )
						{
							print REDFASTA "\n", $line, "\n";
							$save = 1;
							$count ++;
							$unique_alleles{$allelenumber} = 0; # set frequency to 0 as check that was copied

						}
					

					}
					else 
					{ $save = 0; } # knows to skip the sequences lines that follow unwanted identifiers
				}
				elsif ( $line =~ /^[A-Z]/ ) # if is a sequence line, copy only is "save" is turned on by wanted identifier
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

-intab: Table as csv with rows as loci.
-din: Directory with FASTA files where "locusname.FAS" is the file name, like BACT000001.FAS
-dout: Directory where filtered FASTA files will be saved.
-dup: Can give 
EOU

print "Quit because: $message\n";
print "ARGV was " . join (", ", @ARGV) . "\n";
}
