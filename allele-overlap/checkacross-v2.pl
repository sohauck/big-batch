#!/usr/bin/perl

#--------------------------------------------------------
# PROGRAM: checkacross-v2.pl
# AUTHOR:  Sofia Hauck
# CREATED: 20.06.2014
# UPDATED: ----------
# VERSION: v1.00
#--------------------------------------------------------
# VERSION HISTORY
# v1.00 (20.06.2014) created from spilicing older perl script together
# v2.00 (09.04.2014) major re-write, again splicing things 
#--------------------------------------------------------
# transposing adapted from http://stackoverflow.com/questions/3249508/transpose-in-perl
#--------------------------------------------------------

use strict;
use warnings;

# Declares subroutines
sub Usage( ; $ );

# Defines scalars needed from command line
my $fIn1;
my $fIn2;
my $fOut;
my $skipZero;

# Get Command line options, exits if conditions don't look right
if( scalar(@ARGV) < 1 ) { Usage("Not enough command line options"); exit; }
my $i = 0;
my $arg_cnt = 0; 
for ($i=0; $i<=$#ARGV; $i++)
{
	if($ARGV[$i] eq "-h")	       { Usage("You asked for help"); exit; }
	if($ARGV[$i] eq "-in1")         { $fIn1 = $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-in2")         { $fIn2 = $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-out")         { $fOut = $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-skip")         { $skipZero = $ARGV[$i+1] || ''; $arg_cnt++; }
}


# Command line option checks & file checks
if(! defined $fIn1) { Usage("Missing Option: -in1 <FILE>"); exit; }
if(! defined $fIn2) { Usage("Missing Option: -in2 <FILE>"); exit; }
if(! defined $fOut) { Usage("Missing Option: -out <FILE>"); exit; }

if(! -e $fIn1)   	{ Usage("First input file does not exist: $fIn1"); exit; }
if(! -e $fIn2)   	{ Usage("Second input file does not exist: $fIn2"); exit; }
if(-e $fOut) 		{ Usage("Outfile file already exists: $fOut"); exit; }



#creates empty hash table for file information to go
# goes through each file, puts their information into %locitable
# key is locus name, value is an array with all alleles, so multiple appearances mean is shared
my %locitable = ();
readInfile ($fIn1);
readInfile ($fIn2);

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



sub readInfile
{
	my $infile = $_[0];
	open(INFILE, $infile) or die "Cannot open $infile\n";
		while ( my $line = <INFILE> )											
		{        
			# prepare line by removing line break and making sure is in csv format
			chomp($line); 		
			$line =~ s/\r\n?|\n//g;
      			$line =~ s/\t/,/g;
      			
      			# each allele number into an element of array, locus name at $locus
			my @linedata = split(/,/, $line);
			my $locus = shift(@linedata);
			
			# adds every allele number to loci table under that locus name, except "ignore" values
			foreach my $allele (@linedata)
			{
				if (length($skipZero) > 0 && $allele == "0")
				{ } #skipping it since ignored value 
				else
				{ push @{$locitable{$locus}}, $allele; }
			}	
		}  
	close(INFILE);
}

#---------------------------------------------------------------
# Subroutines
#---------------------------------------------------------------

sub Usage( ; $ )
{
	my $message = $ARGV[0] || '';

	print << 'EOU';

Description:
  Checks if alleles are shared over two groups

Example Input Format:
  BACT01	1,23	
  BACT02	0,12,2	

Example Output Format:
  BACT01	specific	2	2	0	(blank)	
  
Usage:

-in1    <FILE> - input filename
-in2	<FILE> - input filename
-out	<FILE> - output filename
-ign	Value to ignore, usually 0 for missing locus
-h     - print usage instructions

EOU

print "Quit because: $message\n";
print "ARGV was " . join (", ", @ARGV);
}