#!/usr/bin/perl

#--------------------------------------------------------
# PROGRAM: gene-diversity.pl
# AUTHOR:  Sofia Hauck
# CREATED: 05.02.2016
#--------------------------------------------------------


use strict;
use warnings;
use LWP::Simple;

$| = 1; # for dynamic output

# Declares subroutines
sub Usage( ; $ );

# Get Command line options, exits if conditions don't look right
my $fTable;			# table with isolates and loci
my $transpose; 
my $dFAS;			# directory where all the FASTA files for each locus are, if have yet
my $dup = 1;		# turned to the smallest number of duplicates necessary for locus to be considered, default is 1



if( scalar(@ARGV) < 2 ) { Usage("Not enough command line options"); exit; }
my $i = 0; 
my $arg_cnt = 0; 
for ($i=0; $i<=$#ARGV; $i++)
{
	if($ARGV[$i] eq "-h")		{ Usage("You asked for help"); exit; }
	if($ARGV[$i] eq "-table")		{ $fIn  = $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-transpose")		{ $dOut = $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-FAS")		{ $fTab = $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-dup")		{ $dup  = $ARGV[$i+1] || ''; $arg_cnt++; }
}

# Command line option checks
if(! defined $fIn) { Usage("Missing Option: input file <FILE>"); exit; }
if(! -e $fIn) { Usage("Input file does not exist: $fIn"); exit; }
if(! defined $fTab)  { Usage("Missing Option: -table <FILE>"); exit; }
if(! defined $dFAS)  { Usage("Missing Option: -dfasta DIRECTORY"); exit; }
if(! defined $dOut)  { Usage("Missing Option: -out DIRECTORY"); exit; }

if(! -e $fTab)  { Usage("Input file does not exist: $fTab"); exit; }
if(! -e $dFAS)  { Usage("Input directory doesn't exist: $dFAS"); exit; }
if(  -e $dOut)  { Usage("Output directory already exists: $dOut"); exit; }



# Preparing for folder for output
if( -e $dOut) { Usage("Output directory already exists: $dOut"); exit; }
mkdir $dOut; 

print "Taking alleles from BIGSDB via the API, now up to...\n";


# Makes output directory and prepares output count file, including header
mkdir $dOut;
my $uniquefile = $dOut . "-count-nuc.txt";

if( -e $uniquefile)  { Usage("Output file already exists: $uniquefile"); exit; }

open(UNIQUENUC, '>', $uniquefile) or die "Cannot open $uniquefile\n";
	print UNIQUENUC "locus,count-nuc,missing\n";


# Transposing, starting with set up...
print "\nTransposing...";
my @original = (); #where data will go in the beginning
my @transposed = (); #where data is moved as its processed
my $columncount = 0; my $rowcount = 0;

# reads table into the @original array
open(INFILE, $fTab) or die "Cannot open $fTab\n";
	while ( my $line = <INFILE> )		
	{
		chomp($line); $line =~ s/\r\n?|\n//g; #just in case, removes all other types of line break
	
		my @row = split (/\t/, $line); # separates by tabs

		$rowcount = @row;
		$columncount ++;
		push ( @original, [ @row ] );
	} 
close(INFILE);

if ( $columncount == 1 ) { die "You only have one column so your csv is probably not in Unix (LF) format.\n"; }



# Read in the main table

my $tableref = ReadTableIn ( $fTable );
my @aoaTable = @$tableref; #puts table, in array of arrays format, into aoaTable
if ( $transpose eq "Yes" )
{
	$tableref = TransposeTable ( \@aoaTable );
	@aoaTable = @$tableref;
}




# Copying only relevant loci
print "Now filter-copying...\n";

foreach my $locusrow (@aoaTable) # loop per locus
{
	# empty hash to find unique values, set locus name
	my %unique_alleles = ();
	
	my $locusname = shift(@{$locusrow});
		
	# populates that hash with unique allele numbers as keys and their number of appearances as values
	foreach my $allele (@{$locusrow}) #go through each allele in locus
    	{
    		if ( $allele =~ /;/ ) # in case of paralogous loci
    		{  
    			my @paralogous = split (/;/, $allele);
    			    			
    			foreach my $paraallele (@paralogous) #go through each allele in locus
				{ $unique_alleles{$paraallele}++; } #increase the frequency count	
    		}
    		    		
    		else # 
    		{ $unique_alleles{$allele}++; } 
	}
	
	my $originalFAS = $dFAS."/".$locusname.".FAS";
	my $reducedFAS = $dOut."/".$locusname.".FAS";

	if ( open(FULLFASTA, $originalFAS) )
	{
		open(REDFASTA, '>', $reducedFAS) or die "Cannot open $reducedFAS\n";
		{       
			my $save = 1; # whether to copy the sequence following a >identifier line
			my $count = 0; # counts how many alleles are copied over, for count-nuc file
			while ( my $line = <FULLFASTA> )
			{
				if ( $line =~ /^>/ ) # if line is >identifier
				{
					chomp $line; 
					my @alleletitle = split ('_', $line); #if in BIGS BACT01_1 format
					my $allelenumber = pop @alleletitle;
					$allelenumber =~ s/[>]+//g; #in case > is still in there if there is no other ID
					
					if ( exists($unique_alleles{$allelenumber}) ) # exists in hash of wanted loci
					{

						if ($unique_alleles{$allelenumber} >= $dup ) # and in frequence at or above cutoff
						{
							print REDFASTA "\n", $line, "\n";
							$save = 1;
							$count ++;
							$unique_alleles{$allelenumber} = 0; # set frequency to 0 as check that was copied
						}
						else 
						{ $save = 0; } # knows to skip the sequences lines that follow unwanted identifiers
					}
					else { $save = 0; }
				}
				
				elsif ( $line =~ /^[A-Z]/ ) # if is a sequence line, copy only is "save" is turned on by wanted identifier
				{
					chomp $line; 
					if ($save == 1)    { print REDFASTA $line; }
					elsif ($save == 0) { }
				}	
			}
			
			my $missing = 0; # sets the missing count back to empty
			
			if ( exists($unique_alleles{"0"}) ) # move count away from empty if any missing alleles were seen
			{ $missing = $unique_alleles{"0"}; } # gives the value in the frequency hash when the key is allele "0", the missing allele
			print UNIQUENUC $locusname . "," . $count . "," . $missing . "\n"; 

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






close(UNIQUENUC);
print "All done!\n";

#---------------------------------------------------------------
# Subroutines
#---------------------------------------------------------------

sub Unique
 {
    my %seen;
    grep !$seen{$_}++, @_;
}

sub ReadTableIn
{
	my $infile = $_[0];
	my @aoaTable = ();
	my $rowcount = 0; 
	
	open(INTABLE, $infile) or die "Cannot open $infile\n";

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

sub GetFASTASeqs
{
	my $dbname	= $_[0];
	my @loci	= @$_[1];

	my $dOut = #??
	chdir $dOut;
	
	# For each locus in the list, get URL of FASTA file and save it to directory
	foreach my $locus ( @loci )											
	{	
		my $url = "http://rest.pubmlst.org/db/".$dbname."/loci/".$locus."/alleles_fasta";
	
		my $file = $locus.".FAS";
		getstore($url, $file);
	
		print "\r$line";
	}

	print "Now completed importing FASTA sequences. They are at: $dOut\n";

}

sub Translate 
{

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


sub Usage( ; $ )
{
	my $message = $_[0] || '';

	print << 'EOU';

keepseenloci.pl

Description:
	  
Usage:
keepseenloci.pl [ options ]

Copies only the alleles in FASTA files that appear in a locus/isolate table. 

-tin: Table as csv with columns as loci.
-din: Directory with FASTA files where "locusname.FAS" is the file name, like BACT000001.FAS
-dout: Directory where filtered FASTA files will be saved.
-dup: Can set a mininum frequence for allele appearance (default is 1) 
EOU

print "\n\nQuit because: $message\n\n";
print "ARGV was " . join (", ", @ARGV) . "\n";
}
