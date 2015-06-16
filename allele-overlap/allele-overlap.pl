#!/usr/bin/perl

#--------------------------------------------------------
# PROGRAM: allele-overlap.pl
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
my $fIn1;
my $fIn2;
my $fOut;

# Get Command line options, exits if conditions don't look right
if( scalar(@ARGV) < 1 ) { Usage("Not enough command line options"); exit; }
my $i = 0;
my $arg_cnt = 0; 
for ($i=0; $i<=$#ARGV; $i++)
{
	if($ARGV[$i] eq "-h")	       { Usage("You asked for help"); exit; }
	if($ARGV[$i] eq "-in1")         { $fIn = $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-in2")         { $fIn = $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-out")        { $fOut = $ARGV[$i+1] || ''; $arg_cnt++; }
}

# Command line option checks
if(! defined $fIn) { Usage("Missing Option: -i|--in <FILE>"); exit; }
if(! defined $fOut) { Usage("Missing Option: -o|--out <FILE>"); exit; }

# File checks
if(! -e $fIn)   { Usage("Input file does not exist: $fIn"); exit; }
if(  -e $fOut)  { Usage("Output file already exists: $fOut"); exit; }


# Transponsing if necessary
if($transpose == 1)
{
	my $t_fIn = $fIn . "-transp"; 
	my @original = (); #where data will go in the beginning
	my @transposed = (); #where data will go after being transposed
	my $columncount = 0; #refers to the final row count and column count, not the initial
	my $rowcount = 0;
	

	open(INFILE, $fIn) or die "Cannot open $fIn\n";
		while ( my $line = <INFILE> )		
		{
		   chomp($line);
		   $line =~ s/\r\n?|\n//g; #just in case, removes all other types of line break
		   $line =~ s/,/\t/g; #if in comma format, switches it to tabs
	   
		   my @row = split ('\t', $line);
		   if ( $rowcount < 2 ) #get rowcount only from first row, since will have no missing values
			{ $rowcount = @row; }
		   $columncount ++; #add to columncount with each row (will become column)
		   push ( @original, [ @row ] );
		} 
	close(INFILE);

	for my $row (@original) #actually does the transposition
	{
		for my $column (0 .. $rowcount) 
		{ push(@{$transposed[$column]}, $row->[$column]); }
	}
	
	@transposed = splice (@transposed, 0, $rowcount); #removes empty rows if more isolates than loci

	open(TRANSP, '>>', $t_fIn) or die "Cannot open $t_In \n"; #opens your results file as appending to it
		for my $new_row (@transposed) 
		{
			for my $new_col (@{$new_row}) 
			{
				print TRANSP $new_col, "\t";
			}
			print TRANSP "\n";
		}
	close(TRANSP);
}


#from this point on, if asks for infile, check if transpose, then use middle file?

if ($transpose == 1)
{ open(INFILE, $t_fIn) or die "Cannot open $t_fIn \n"; }
else
{ open(INFILE, $fIn) or die "Cannot open $fIn \n"; }

	<INFILE>; # Skips header
	while ( my $line = <INFILE> ) #goes through each row in the file at a time			
	{
		$line =~ s/\r\n?|\n//g; 		#removes line breaks
		$line =~ s/,/\t/g; 			#switches to tab separators
        
		my @linedata = grep {length()} split(/,/, $line); #put each allele into one item of an array, skips empty cells 

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
