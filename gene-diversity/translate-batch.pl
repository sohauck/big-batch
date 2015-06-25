#!/usr/bin/perl

#--------------------------------------------------------
# PROGRAM: translate-batch.pl
# AUTHOR:  Sofia Hauck
# CREATED: 15.06.2015
#--------------------------------------------------------

use strict;
use warnings;
$| = 1;

# Declares subroutines
sub Usage( ; $ );

# Get command line options & file checks, including creating output directory and output count file
my $dir = shift(@ARGV);
if(! defined $dir) { Usage("Missing Option: a directory as the first argument"); exit; }
if(! -e $dir)      { Usage("Input directory does not exist: $dir"); exit; }

my @splitdir =   split(/\//, $dir);
my $oldfolder =  pop @splitdir;
my $transdir =   join("/", @splitdir) . "/" . $oldfolder . "-t/"; 
my $uniquefile = join("/", @splitdir) . "/" . $oldfolder. "-count-aa.txt";

if( -e $transdir)     { Usage("Output directory already exists: $transdir"); exit; }
if( -e $uniquefile)   { Usage("Output file already exists: $uniquefile"); exit; }

mkdir $transdir;
open(UNIQUEAA, '>', $uniquefile) or die "Cannot open $uniquefile\n"; # setting up header of output count file 
	print UNIQUEAA "locus,count-aa,min-length,max-length,avg-length\n";


# Get names of all the files in the directory
opendir (ORIGDIR, $dir) or die "Cannot open directory: $!";
	my @files = readdir ORIGDIR;
	@files = grep(/^([A-Z]|[a-z]|[0-9])/,@files);
	if ($#files < 1)
	{ die "It looks like you have no FASTA files to align. If their file names don't begin with letters or numbers they are not being included.\n";}
closedir ORIGDIR;


# Loop that aligns translates all the files and puts them in new directory
print "Translated up to...\n";

foreach my $file (@files)
{
 	my $infile = $dir . "/" . $file; # original file from folder given in command line
 	my $outfile = $transdir . $file; # aligned file into translated folder, with same file name
 	
	my %unique = (); # empty hash to check for synonymous mutations 

	open(NUCLEOTIDE, $infile) or die "Cannot open $infile\n";
	open(AMINOACID, '>', $outfile) or die "Cannot open $outfile\n";


 	while ( my $line = <NUCLEOTIDE> )
	{
		if ( $line =~ /^>/ ) # if line is title of sequence
		{ 
			print AMINOACID $line;
		}
		elsif ( $line =~ /^[A-Z]/ ) # if line is sequence
		{
			chomp $line;
			my $peptide = translate($line, "1");
			
			# for later checking the count of unique amino acids
	    		if ( !exists($unique{$peptide})) #if allele doesn't exist in unique-hash
    			{ $unique{$peptide} = 1; } #then add it with allele number as key, value as 1
			
			print AMINOACID $peptide, "\n";
		}	
	}
	
 	my @filearray = split (/\./, $file);
 	my $locusname = shift @filearray;
 	print "\r$locusname";
 
 	# counting number of unique amino acids
 	my $count;
 	my @lengths;
 	
 	foreach my $key ( sort keys %unique ) 
 	{ 	$count++;
 		push (@lengths, length($key));
 	}
 	
 	use List::Util qw( min max sum );
	my $min = min @lengths;
	my $max = max @lengths;
	my $avg = sum @lengths / ($#lengths + 1);
 	
 	print UNIQUEAA  $locusname .",". $count .",". $min .",". $max .",". $avg ."\n";
 	 
 	
 	close NUCLEOTIDE; close AMINOACID;
 
 	
 	# So you have something to watch while it runs... 
 	# print "\r$locusname";
}

print "\nComplete! Results are at $transdir.\n";

close UNIQUEAA;


#---------------------------------------------------------------
# Subroutines
#---------------------------------------------------------------

sub Usage( ; $ )
{
	my $message = $_[0] || '';

	print << 'EOU';

translate-batch.pl

Description:
	  
Usage:
translate-batch.pl [ options ]

Takes FASTA files from a folder, translates them into amino acid sequences.
Also reports number of unique amino acid sequences (useful as comparing to unique nucleotide sequences).

First argument is folder with FASTA files, everything else gets ignore.

EOU

print "Quit because: $message\n";
print "ARGV was " . join (", ", @ARGV) . "\n";
}

sub translate  {

# Author: John Nash
#
# Copyright (c) Government of Canada, 2000-2012,
#   all rights reserved.
#
# Licence: This script may be used freely as long as no fee is charged 
#   for use, and as long as the author/copyright attributions 
#   are not removed.


# Does:  Translates incoming string from nucleic acid to peptide
# Expects: nucleic acid string ($sequence_str); 
#          reading frame value: +1, +2, +3, -1, -2, -3
# Returns: string of translated peptide
# Uses:
    
# vars for translation:
    my ($codon, $peptide, $seq_str, $frame);
    
# reading frame:
    my $startpoint;
    
# incoming:
    $seq_str = uc $_[0];
    
# assess the reading frame:
    $frame = $_[1];
    $seq_str =~ tr/natgcbdkrvhmyxwsNATGCBDKRVHMYXWS
    /ntacgvhmybdkrxswNTACGVHMYBDKRXSW/ if ($frame < 0);
    $startpoint = (abs $frame) - 1;
    
# Perform the translation:
    for (my $j=$startpoint; $j < length $seq_str; $j+=3) {
	$codon=(substr($seq_str, $j, 3));
	if (($codon eq "GCA") || 
	    ($codon eq "GCC") ||
	    ($codon eq "GCG") ||
	    ($codon eq "GCT") ||
	    ($codon eq "GCR") ||
	    ($codon eq "GCY") ||
	    ($codon eq "GCN")) {
	    $peptide.="A";
	}   # Ala A
	
	elsif (($codon eq "CGA") || 
	       ($codon eq "CGC") ||
	       ($codon eq "CGG") ||
	       ($codon eq "CGT") ||
	       ($codon eq "CGR") ||
	       ($codon eq "CGY") ||
	       ($codon eq "CGN") ||
	       ($codon eq "AGA") ||
	       ($codon eq "AGG") ||
	       ($codon eq "AGR")) {
	    $peptide.="R";
	}   # Arg R
	
	elsif (($codon eq "AAC") || 
	       ($codon eq "AAT") ||
	       ($codon eq "AAY")) {
	    $peptide.="N";
	}   # Asn N
	
	elsif (($codon eq "GAC") || 
	       ($codon eq "GAT") ||
	       ($codon eq "GAY")) {
	    $peptide.="D";
	}   # Asp D
	
	elsif (($codon eq "TGC") || 
	       ($codon eq "TGT") ||
	       ($codon eq "TGY")) {
	    $peptide.="C";
	}   # Cys C
	
	elsif (($codon eq "CAA") || 
	       ($codon eq "CAG") ||
	       ($codon eq "CAR")) {
	    $peptide.="Q";
	}   # Gln Q
	
	elsif (($codon eq "GAA") || 
	       ($codon eq "GAG") ||
	       ($codon eq "GAR")) {
	    $peptide.="E";
	}   # Glu E
	
	elsif (($codon eq "GGA") || 
	       ($codon eq "GGC") ||
	       ($codon eq "GGG") ||
	       ($codon eq "GGT") ||
	       ($codon eq "GGR") ||
	       ($codon eq "GGY") ||
	       ($codon eq "GGN")) {
	    $peptide.="G";
	}   # Gly G
	
    elsif (($codon eq "CAC") || 
	   ($codon eq "CAT") ||
	   ($codon eq "CAT")) {
	$peptide.="H";
    }   # His H
	
	elsif (($codon eq "ATA") || 
	       ($codon eq "ATC") ||
	       ($codon eq "ATT") ||
	       ($codon eq "ATH")) {
	    $peptide.="I";
	}   # Ile I
	
	elsif (($codon eq "CTA") || 
	       ($codon eq "CTC") ||
	       ($codon eq "CTG") ||
	       ($codon eq "CTT") ||
	       ($codon eq "CTR") ||
	       ($codon eq "CTY") ||
	       ($codon eq "CTN") ||
	       ($codon eq "TTA") ||
	       ($codon eq "TTG") ||
	       ($codon eq "TTR")) {
	    $peptide.="L";
	}   # Leu L
	
	elsif (($codon eq "AAA") || 
	       ($codon eq "AAG") ||
	       ($codon eq "AAR")) {
	    $peptide.="K";
	}   # Lys K
	
	elsif (($codon eq "ATG")) {
	    $peptide.="M";
	}   # Met M
	
	elsif (($codon eq "TTC") || 
	       ($codon eq "TTT") ||
	       ($codon eq "TTY")) {
	    $peptide.="F";
	}   # Phe F
	
	elsif (($codon eq "CCA") || 
	       ($codon eq "CCC") ||
	       ($codon eq "CCG") ||
	       ($codon eq "CCT") ||
	       ($codon eq "CCR") ||
	       ($codon eq "CCY") ||
	       ($codon eq "CCN")) {
	    $peptide.="P";
	}   # Pro P
	
	elsif (($codon eq "TCA") || 
	       ($codon eq "TCC") ||
	       ($codon eq "TCG") ||
	       ($codon eq "TCT") ||
	       ($codon eq "TCR") ||
	       ($codon eq "TCY") ||
	       ($codon eq "TCN") ||
	       ($codon eq "AGC") ||
	       ($codon eq "AGT") ||
	       ($codon eq "AGY")) {
	    $peptide.="S";
	}   # Ser S
	
	elsif (($codon eq "ACA") || 
	       ($codon eq "ACC") ||
	       ($codon eq "ACG") ||
	       ($codon eq "ACT") ||
	       ($codon eq "ACR") ||
	       ($codon eq "ACY") ||
	       ($codon eq "ACN")) {
	    $peptide.="T";
	}   # Thr T
	
    elsif (($codon eq "TGG")) {
	$peptide.="W";
    }   # Trp W
	
	elsif (($codon eq "TAC") || 
	       ($codon eq "TAT") ||
	       ($codon eq "TAY")) {
	    $peptide.="Y";
	}   # Tyr Y
	
	elsif (($codon eq "GTA") || 
	       ($codon eq "GTC") ||
	       ($codon eq "GTG") ||
	       ($codon eq "GTT") ||
	       ($codon eq "GTR") ||
	       ($codon eq "GTY") ||
	       ($codon eq "GTN")) {
	    $peptide.="V";
	}   # Val V
	  
	elsif (($codon eq "TAA") || 
	       ($codon eq "TAG") ||
	       ($codon eq "TAR") ||
	       ($codon eq "TGA")) {
	    $peptide.="*";
	}   # Stop *
	
	else {
	    ;
	} # do nothing for now...
    } # end of  for ($j=0; etc...
    
# return the value:
    return $peptide;
}