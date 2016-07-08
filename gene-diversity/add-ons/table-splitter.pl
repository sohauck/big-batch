#!/usr/bin/perl

#--------------------------------------------------------
# PROGRAM: table-splitter.pl
# AUTHOR:  Sofia Hauck
# CREATED: 18.05.2016
#--------------------------------------------------------

use strict;
use warnings;
use List::Util qw( min max sum ); # for measuring lengths

$| = 1; # for dynamic output
my $version = "0.9";

# Declares subroutines
sub Usage( ; $ );



# Get Command line options, exits if conditions don't look right
my $fTable;			# file where first line is a header, others are individual records, ideally randomised already
my $method = "";			# whether table will be partitioned and broken into accumulation steps
my $seccount = 0;		# how many files the main one will eventually be split into
my $dOut; 			# directory where results will go 

my $i = 0; 
my $arg_cnt = 0; 
for ($i=0; $i<=$#ARGV; $i++)
{
 	if($ARGV[$i] eq "-help")		{ Usage("You asked for help"); exit; }
 	if($ARGV[$i] eq "-table")		{ $fTable		= $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-partition")	{ $method = "partition"; $seccount = $ARGV[$i+1]; $arg_cnt++; }
	if($ARGV[$i] eq "-accumulate")	{ $method = "accumulate"; $seccount = $ARGV[$i+1]; $arg_cnt++; }
	if($ARGV[$i] eq "-dout")		{ $dOut = $ARGV[$i+1]; $arg_cnt++; }
}


# a table of loci and isolates 
if(! defined $fTable) { Usage("Input table file not defined"); exit; }
if(! -e $fTable)  { Usage("Input table file does not exist: $fTable"); exit; }


# Check options make sense
if ( $method eq "" || $seccount == "0" )
{ Usage("You didn't choose a method and an integer for the number of parts: $method and $seccount"); exit;  }


# a directory where results can go
if(! defined $dOut)  { 	Usage("Output directory not defined"); exit; }
if(  -e 	 $dOut)  {  Usage("Output directory already exists: $dOut"); exit; } 

if ( mkdir $dOut )
{ print "\nGetting started... "; }
else { Usage("Could not create output directory: $dOut"); exit; }


# Actually starting the work now


# Reading the table in
my @aTable = ();
	
open(INTABLE, $fTable) or die "Cannot open $fTable\n";

while ( my $line = <INTABLE> )
{
	chomp($line); $line =~ s/\r\n?|\n//g; # clean up line breaks
	push ( @aTable, $line );
} 

my $rowcount = scalar (@aTable); 
if ( $rowcount <= 1 ) { die "You only have one column so something has probably gone wrong with line breaks.\n"; }

close(INTABLE);
	



# pop the header onto its own 
my $headerrow = shift @aTable;
$rowcount --; 

# calculate the sizes of the sections
my $rowspersection = $rowcount / $seccount;
my $remainder = $rowcount % $seccount; 


print "Header begins " . substr ( $headerrow , 0, 20 ). "\n"; 
print "Number of rows: $rowcount \n"; 
print "Sections to split into: is $seccount \n\n"; 

print "rowspersection is $rowspersection \n"; 
print "remainder is $remainder \n\n\n"; 



my @accTable; #empty table to push the already included rows into for accumulated sections

for ( my $section = 1; $section <= $seccount; $section++ ) # loop by the number of sections
{
	my $sectionfile = $dOut . "/" . $method . "-" . $section . ".txt";
			
	open (SECTION, '>', $sectionfile) or die "Cannot open $sectionfile\n"; # then also open where the reduced one will go
	print SECTION $headerrow . "\n"; 
	
	# if on accumulation method, put all the rows so far in there to begin with
	if ( $method eq "accumulate" ) 
	{ 
		foreach my $accrow ( @accTable )
		{ print SECTION $accrow . "\n"; }  
	}

	# in every case, go through the remaining rows, put as many in there are takes in each section
	for ( my $rowsadded = 1; $rowsadded < $rowspersection; $rowsadded++ )
	{
		my $row = shift @aTable;
		print SECTION $row . "\n"; 
		
		if ( $method eq "accumulate" )
		{ push ( @accTable, $row ); }
	}
	
	if ( $section <= $remainder )
	{
		my $row = shift @aTable;
		print SECTION $row . "\n"; 
		
		if ( $method eq "accumulate" )
		{ push ( @accTable, $row ); }
	}
	
	close SECTION; 
}



#---------------------------------------------------------------
# Subroutines
#---------------------------------------------------------------

sub Usage( ; $ )
{
	my $message = $_[0] || '';

	print << 'EOU';

title

Description:  
	  
Usage:
file name

Command line options:
-table: Table as tab-separated, with loci and isolates as rows and columns, or reverse if "transpose" option used.
-transpose: "Yes" if the table has loci as columns and isolates as rows, "No" otherwise. 
-partition: followed by number, how many (roughly) pieces to break apart into
-accumulate: follow by number, how many incremental pieces to break apart into
-dout: output directory

EOU

print "\n\nQuit because: $message\n\n";
print "ARGV was " . join (", ", @ARGV) . "\n";
}
