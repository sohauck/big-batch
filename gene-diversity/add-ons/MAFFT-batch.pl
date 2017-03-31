#!/usr/bin/perl

#--------------------------------------------------------
# PROGRAM: MAFFT-batch.pl
# AUTHOR:  Sofia Hauck
# CREATED: 08.07.2016
#--------------------------------------------------------

use strict;
use warnings;


# Declares subroutines
sub Usage( ; $ );



# Get Command line options, exits if conditions don't look right
my $dIn;									# folder where inputs are (GbGDiv results)
my $dOut;									# folder where outputs are going 
my @mafftarg = ("--quiet");	# start with at least these, can ask for more, add "-mafft --auto" to keep silent

my $i = 0; 
my $arg_cnt = 0; 
for ($i=0; $i<=$#ARGV; $i++)
{
 	if($ARGV[$i] eq "-help")		{ Usage("You asked for help"); exit; }
 	if($ARGV[$i] eq "-dIn")			{ $dIn	= $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-dOut")		{ $dOut = $ARGV[$i+1] || ''; $arg_cnt++; }
 	if($ARGV[$i] eq "-mafft")		{ push (@mafftarg, $ARGV[$i+1]) ; $arg_cnt++; }
}


# a directory where results can go
if(! defined $dIn)  { 	Usage("Input directory not defined"); exit; }
if(! -e 	 $dIn)  {  Usage("Input directory doesn't exist: $dIn"); exit; } 

if(! defined $dOut)  { 	Usage("Output directory not defined"); exit; }
if(  -e 	 $dOut)  {  Usage("Output directory already exists: $dOut"); exit; } 

if ( mkdir $dOut )
{ print "\nGetting started... "; }
else { Usage("Could not create output directory: $dOut"); exit; }


# read all names from files


opendir (OBSDIR, $dIn ) or die "Cannot open directory: $!";
	my @files = readdir OBSDIR;
	@files = grep(/^([A-Za-z0-9])/,@files); # remove any hidden files
	if ($#files < 1)
	{ die "Couldn't find any FASTA files. Maybe their names don't begin with a letter or number?\n";}
	@files = sort @files; # so they're in order and the "Currently up to" is a bit more informative
closedir OBSDIR;



my $preMA;
my $postMA;

for my $file ( @files )
{
 	# where all my files at
 	my $preMA	=	$dIn  ."/". $file; 
 	my $postMA	=	$dOut ."/". $file;

	system ( "mafft " . join (" ", @mafftarg) . " " . $preMA . " > " . $postMA  ); 
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
