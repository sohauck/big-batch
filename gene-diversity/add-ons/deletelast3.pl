#!/usr/bin/perl

#--------------------------------------------------------
# PROGRAM: deletelast3.pl
# AUTHOR:  Sofia Hauck
# CREATED: 05.02.2016
#--------------------------------------------------------

use strict;
use warnings;


# Declares subroutines
sub Usage( ; $ );



# Get Command line options, exits if conditions don't look right
my $dIn;			# folder where inputs are (GbGDiv results)
my $dOut;			# folder where outputs are going 

my $i = 0; 
my $arg_cnt = 0; 
for ($i=0; $i<=$#ARGV; $i++)
{
 	if($ARGV[$i] eq "-help")		{ Usage("You asked for help"); exit; }
 	if($ARGV[$i] eq "-in")			{ $dIn	= $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-out")			{ $dOut = $ARGV[$i+1] || ''; $arg_cnt++; }
}


# a directory where results can go
if(! defined $dIn)  { 	Usage("Input directory not defined"); exit; }
if(! -e 	 $dIn)  {  Usage("Input directory doesn't exist: $dIn"); exit; } 

if(! defined $dOut)  { 	Usage("Output directory not defined"); exit; }
if(  -e 	 $dOut)  {  Usage("Output directory already exists: $dOut"); exit; } 

if ( mkdir $dOut )
{ print "\nGetting started... "; }
else { Usage("Could not create output directory: $dOut"); exit; }


# opening the directory where the observed alleles are and getting all the file names
opendir ( OBSDIR, $dIn ) or die "Cannot open directory: $!";
	my @files = readdir OBSDIR;
	@files = grep(/^([A-Za-z0-9])/,@files); # remove any hidden files
	if ($#files < 1)
	{ die "Couldn't find any FASTA files. Maybe their names don't begin with a letter or number?\n";}
	@files = sort @files; # so they're in order and the "Currently up to" is a bit more informative
closedir OBSDIR;


for my $file ( @files )
{
	my $original = $dIn . "/" . $file;
	my $chopped  = $dOut . "/" . $file;
	
	open(ORIGINAL, $original) ;
		open(CHOPPED, '>', $chopped) or die "Cannot open $chopped\n"; # then also open where the reduced one will go
		
		while ( my $line = <ORIGINAL> ) # reading through original FASTA
		{
			if ( $line =~ /^>/ ) # if line is >identifier
			{	print CHOPPED $line ;	}
		
			elsif ( $line =~ /^([A-Za-z*-])/ ) # if is a sequence line (starts with a letter, - or *) 
			{
				chomp $line; 
				
				$line = substr ( $line, 0, -3 ); 
				
				print CHOPPED $line . "\n" ;
			}	
		}

	close(CHOPPED);
	close(ORIGINAL);
	
}

print "\nDone!!\n";

#---------------------------------------------------------------
# Subroutines
#---------------------------------------------------------------


sub Usage( ; $ )
{
	my $message = $_[0] || '';

	print << 'EOU';

-in
-out 

EOU

print "\n\nQuit because: $message\n\n";
print "ARGV was " . join (", ", @ARGV) . "\n";
}
