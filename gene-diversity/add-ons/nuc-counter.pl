#!/usr/bin/perl

#--------------------------------------------------------
# PROGRAM: nuc-counter.pl
# AUTHOR:  Sofia Hauck
# CREATED: 19.05.2016
#--------------------------------------------------------

use strict;
use warnings;

$| = 1; # for dynamic output

# Declares subroutines
sub Usage( ; $ );


# Get Command line options, exits if conditions don't look right
my $dTable;			# directory with tables of isolates and loci
my $transpose;		# whether the table needs to be transposed or no
my $dOut;			# directory where the results will go
my $dup = 1;		# the smallest number of duplicates necessary for locus to be considered "seen", default is 1


my $i = 0; 
my $arg_cnt = 0; 
for ($i=0; $i<=$#ARGV; $i++)
{
 	if($ARGV[$i] eq "-help")		{ Usage("You asked for help"); exit; }
 	if($ARGV[$i] eq "-table")		{ $dTable		= $ARGV[$i+1] || ''; $arg_cnt++; }
 	if($ARGV[$i] eq "-transpose")	{ $transpose 	= $ARGV[$i+1] || ''; $arg_cnt++; }
 	if($ARGV[$i] eq "-out")			{ $dOut  		= $ARGV[$i+1] || ''; $arg_cnt++; }
 	if($ARGV[$i] eq "-dup")			{ $dup  		= $ARGV[$i+1] || ''; $arg_cnt++; }
}


# Do we have everything we need to start?

# a table of loci and isolates 
if(! defined $dTable) 
{ 
	print "Where is the directory with tab-separated tables of isolates and loci?\n"; 
	$dTable = <STDIN>; 
	chomp $dTable; $dTable =~ s/\s+$//; # removes white spaces and line breaks
}
if(! -e $dTable)  { Usage("Input table directory does not exist: $dTable"); exit; }


# and whether to transpose that table or no
if (! defined $transpose )
{
	print "\n\nDoes the table have loci as columns? \nIf yes, please enter 'Yes' as the table needs transposing. " .
	"\nIf no, that is, if the each row of your table is a locus, transposing is not needed, so enter 'No' below.\n"; 
	$transpose = <STDIN>; 
	if ( ! $transpose =~ /^[YyNn]/ ) # as long as it starts with a "y" or "n"
	{ 
		print "That didn't look like a yes or a no... Try again? \n";
		$transpose = <STDIN>; 
	}
	elsif ( $transpose =~ /^[Yy]/ )
	{ $transpose = "Yes" }
	elsif ( $transpose =~ /^[Nn]/ )
	{ $transpose = "No" }
	else 
	{ Usage("Something went wrong with the transposing options"); exit; }
}


# a directory where results can go
if(! defined $dOut)  
{ 	
	print 	"\n\nName the directory where you'd like the results to go.\n";

	$dOut = <STDIN>; 
	chomp $dOut; $dOut =~ s/\s+$//;

	if(  -e $dOut)  { Usage("Output directory already exists: $dOut"); exit; } 
}



# Now that we've asked for everything, let's check that it's all what we expect... 
print "\n\n\nLet's check if everything is correct before we start...\n\n";

if ( $transpose =~ /^[Nn]/ )
{ print	"Your tables of isolates as columns and loci as rows are held at $dTable.\n"  }
elsif ( $transpose =~ /^[Yy]/ )
{ print	"Your tables of loci as columns and isolates as rows are held at $dTable.\n"  }



# Create output directory
if( -e $dOut) { Usage("Output directory already exists: $dOut"); exit; }

if ( mkdir $dOut )
{ print "Your results will be in the new directory $dOut\n\n\n"; }
else { Usage("Could not create output directory: $dOut"); exit; }




# Actually starting the work now

opendir (TABDIR, $dTable) or die "Cannot open directory: $!";
	my @tabfiles = readdir TABDIR; 
	@tabfiles = grep(/^([A-Za-z0-9])/,@tabfiles); # remove any hidden files (doesn't start with letter or number)
	if ($#tabfiles < 1)
	{ die "Couldn't find any table files. Maybe their names don't begin with a letter or number?\n";}
	@tabfiles = sort @tabfiles; # so they're in order and the "Currently up to" is a bit more informative
closedir TABDIR;


foreach my $tabfile ( @tabfiles )
{
	# Reading the table in
	my $tableaddress = $dTable . "/" . $tabfile; 
	my $tableref = ReadTableIn ( $tableaddress );
	my @aoaTable = @$tableref; #puts table, in array of arrays format, into aoaTable

	if ( $transpose =~ /^[Yy]/ ) # transposes the table and puts it back in aoaTable
	{
		print "Transposing your table...";
		$tableref = TransposeTable ( \@aoaTable );
		@aoaTable = @$tableref;
		print " and done.\n";
	}

	my $headerrow = $aoaTable[0];
	my $isolatecount = scalar (@$headerrow) - 1; # since the first column is the header

	# Printing the header of the results file
	open ( COUNTNUC, '>', $dOut."/Res-CountNuc-".$tabfile.".txt" ) or die "$dOut /ResultsTable.txt"; # then also open where the reduced one will go

		print COUNTNUC	"GbGDiv simplified results table from '". $tabfile ."' records table.\n" . 
						"Isolate records read from table: ". $isolatecount . "\n\n".
						join ("\t", ("Locus","Missing","Paralogous","CountNuc") ) . "\n" ;


	# Going through table and keeping the observed alleles, plus counting missing and unique nucleotide sequences 
	print "Now reading your table in...\n\n";


	foreach my $locusrow (@aoaTable) # loop per locus
	{
		# empty hash to find unique values, values to be found
		my %unique_alleles = ();
		my $paralogous = 0; 
		my $missing = 0;

		# set locus name
		my $locusname = shift(@{$locusrow});
		$locusname =~ /^(\S+)/; $locusname = $1; # if there are aliases or anything else after the locus ID
			
		# populates that hash with unique allele numbers as keys and their number of appearances as values
		foreach my $allele (@{$locusrow}) #go through each allele in locus
		{
				if ( length ( $allele ) == 0 || $allele eq "X" || $allele eq "I" )
				{	$unique_alleles{"0"}++;	} # if cell is empty or has "X", file as "missing locus"
			
				else
				{
					my @alleles = split (/;/, $allele); # in case of paralogous loci
					if ( scalar(@alleles) >= 2 )
					{	$paralogous++; }

					foreach (@alleles) # go through each allele in locus
					{	$unique_alleles{$_}++;	} #increase the frequency count	
				}
		}

		my $countnuc = scalar ( keys %unique_alleles );

		if ( exists($unique_alleles{"0"}) )
		{	
			$missing = $unique_alleles{"0"}; # gives the value in the frequency hash when the key is allele "0", the missing allele
			$countnuc -- ; # discount one "allele" from count since 0 is not an allele
		} 		

		# find the file where the original FASTA sequences are	
	
		print COUNTNUC join ("\t", ($locusname, $missing, $paralogous, $countnuc) ), "\t\n";	
					
	} # close per locus loop

	close COUNTNUC; # close results file to adding one line per locus as it goes through table-reading loop
	print "\nTable reading complete.\n";

}



#---------------------------------------------------------------
# Subroutines
#---------------------------------------------------------------


sub ReadTableIn
{
	my $infile = $_[0];
	my @aoaTable = ();
	my $rowcount = 0; 
	
	open(INTABLE, $infile) or die "Cannot open the table $infile\n";

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
	my $aoaRef = $_[0]; 
	my @original = @$aoaRef; # @original is now array of references to arrays
		
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
		
	@transposed = splice (@transposed, 0, scalar(@{$original[0]})); #removes empty rows if more isolates than loci
	
	return \@transposed; 
}


sub Usage( ; $ )
{
	my $message = $_[0] || '';

	print << 'EOU';

...

Description: ... 
	  
Usage:
nuc-counter.pl

Command line options:
-table: Table as tab-separated, with loci and isolates as rows and columns, or reverse if "transpose" option used.
-transpose: "Yes" if the table has loci as columns and isolates as rows, "No" otherwise. 
-out: Directory where all results will be saved.
-dup: Value for the mininum frequency for allele appearance to be included. Default is "1".
EOU

print "\n\nQuit because: $message\n\n";
print "ARGV was " . join (", ", @ARGV) . "\n";
}
