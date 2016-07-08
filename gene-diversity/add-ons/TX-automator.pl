#!/usr/bin/perl

#--------------------------------------------------------
# PROGRAM: TX-automator.pl
# AUTHOR:  Sofia Hauck
# CREATED: 08.07.2016
#--------------------------------------------------------

use strict;
use warnings;


# Declares subroutines
sub Usage( ; $ );



# Get Command line options, exits if conditions don't look right
my $chopped;			# folder with raw FASTA without end codon
my $aligned;			# folder with aligned FASTA
my $dOut;				# folder where outputs are going 
my $translatorxloc; 	# where the TranslatorX file is

my $i = 0; 
my $arg_cnt = 0; 
for ($i=0; $i<=$#ARGV; $i++)
{
 	if($ARGV[$i] eq "-help")		{ Usage("You asked for help"); exit; }
 	if($ARGV[$i] eq "-chopped")		{ $chopped	= $ARGV[$i+1] || ''; $arg_cnt++; }
 	if($ARGV[$i] eq "-aligned")		{ $aligned	= $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-out")		{ $dOut = $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-TX")			{ $translatorxloc = $ARGV[$i+1] || ''; $arg_cnt++; }
}


# a directory where results can go
if(! defined $chopped)  { 	Usage("Input directory not defined"); exit; }
if(! -e 	 $chopped)  {  Usage("Input directory doesn't exist: $chopped"); exit; } 

if(! defined $aligned)  { 	Usage("Input directory not defined"); exit; }
if(! -e 	 $aligned)  {  Usage("Input directory doesn't exist: $aligned"); exit; } 


if(! defined $dOut)  { 	Usage("Output directory not defined"); exit; }
if(  -e 	 $dOut)  {  Usage("Output directory already exists: $dOut"); exit; } 

if ( mkdir $dOut )
{ print "\nGetting started... "; }
else { Usage("Could not create output directory: $dOut"); exit; }


# read all names from files


opendir ( OBSDIR, $chopped ) or die "Cannot open directory: $!";
	my @files = readdir OBSDIR;
	@files = grep(/^([A-Za-z0-9])/,@files); # remove any hidden files
	if ($#files < 1)
	{ die "Couldn't find any FASTA files. Maybe their names don't begin with a letter or number?\n";}
	@files = sort @files; # so they're in order and the "Currently up to" is a bit more informative
closedir OBSDIR;



my $nucfile;
my $aafile;
my $TXfile;

for my $file ( @files )
{
	my ($locusname, $extension) = split (/\./, $file);

 	# where all my files at
 	my $nucfile	=	$chopped ."/". $file; 
 	my $aafile	=	$aligned ."/". $file; 
	my $TXfile  =	$dOut."/".$locusname;

	system ( "perl " . $translatorxloc . " -i ". $nucfile ." -a ". $aafile ." -o ". $TXfile );
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
-dIn: folder where inputs are (GbGDiv results)
-dOut: folder where outputs are going 
-TX: where the TranslatorX file is

EOU

print "\n\nQuit because: $message\n\n";
print "ARGV was " . join (", ", @ARGV) . "\n";
}
