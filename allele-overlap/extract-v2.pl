#!/usr/bin/perl

#--------------------------------------------------------
# PROGRAM: extract-v2.pl
# AUTHOR:  Sofia Hauck
# CREATED: 20.06.2014
# UPDATED: ----------
# VERSION: v1.00
#--------------------------------------------------------
# VERSION HISTORY
# v1.00 (20.06.2014) created from spilicing older perl script together
# v2.00 (09.04.2014) major re-write, against splicing things 
#--------------------------------------------------------
# transposing adapted from http://stackoverflow.com/questions/3249508/transpose-in-perl
#--------------------------------------------------------

use strict;
use warnings;

# Declares subroutines
sub Usage( ; $ );

# Defines scalars needed from command line
my $fIn;

# Get Command line options, exits if conditions don't look right
if( scalar(@ARGV) < 1 ) { Usage("Not enough command line options"); exit; }
my $i = 0;
my $arg_cnt = 0; 
for ($i=0; $i<=$#ARGV; $i++)
{
	if($ARGV[$i] eq "-h")	       { Usage("You asked for help"); exit; }
	if($ARGV[$i] eq "-in")         { $fIn = $ARGV[$i+1] || ''; $arg_cnt++; }
}

my $fOut = $fIn . "-unique"; #names outfile


# Command line option checks & file checks
if(! defined $fIn) { Usage("Missing Option: -i|--in <FILE>"); exit; }
if(! -e $fIn)   { Usage("Input file does not exist: $fIn"); exit; }
if(-e $fOut) { Usage("Outfile file already exists: $fOut"); exit; }


open(INFILE, $fIn) or die "Cannot open $fIn\n";
	<INFILE>; #ignores first row (header)
	while ( my $line = <INFILE> ) 
	{
		#prepare line: removes all types of line breaks, switches to csv if tab-separate
		chomp($line);
		$line =~ s/\r\n?|\n//g;
        	$line =~ s/\t/,/g;
        
        	# each allele as element in an array, name of locus in $locus
     		my @linedata = grep {length()} split(/,/, $line);
       		my $locus = shift(@linedata);
		
		# creates hash to make list of unique alleles
    		my %unique_alleles = (); 
    	
		foreach my $allele (@linedata) #goes through each allele in array
		{
			if ( !exists($unique_alleles{$allele})) #if it's not in the hash table
				{ $unique_alleles{$allele} = 1; } #creates it! 
		}
		
		#empties unique array, then add all unique keys from hash into array
		my @unique = ();
		
		foreach my $key ( sort keys %unique_alleles )
		{
		     push (@unique, $key)
		}
						
		# appends results file with name of locus, then entire array of unique alleles
		open(OUTFILE, '>>', $fOut) or die "Cannot open $fOut\n";
		     print OUTFILE $locus . "," . join (',',@unique) . "\n";
		close(OUTFILE);
	} 
close(INFILE);

print 	"All done!\nCheck your results at $fOut. Each line begins with the locus name and contains all unique alleles at that locus.\n";



#---------------------------------------------------------------
# Subroutines
#---------------------------------------------------------------

sub Usage( ; $ )
{
	my $message = $ARGV[0] || '';

	print << 'EOU';

Description:
  Takes two lists, removes all items in second list from first list, returns remaining

Example Input Format:
    		123	124	125	127
  BACT01	1	2	1	1

Example Output Format:
  BACT01	1,2	
  
Usage:

-in    <FILE> - input filename
-h     - print usage instructions

EOU

print "Quit because: $message\n";
print "ARGV was " . join (", ", @ARGV);
}
