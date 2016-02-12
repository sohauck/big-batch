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
my $transpose = "No";	# whether the table needs to be transposed or no
my $dFAS;			# directory where all the FASTA files for each locus are, if they don't need to be exported
my $dOut;			# directory where the results will go
my $dup = 1;		# the smallest number of duplicates necessary for locus to be considered "seen", default is 1


my $i = 0; 
my $arg_cnt = 0; 
for ($i=0; $i<=$#ARGV; $i++)
{
	if($ARGV[$i] eq "-help")		{ Usage("You asked for help"); exit; }
	if($ARGV[$i] eq "-table")		{ $fTable  = $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-transpose")	{ $transpose = "Yes"  || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-FASTA")		{ $dFAS = $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-out")			{ $dOut  = $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-dup")			{ $dup  = $ARGV[$i+1] || ''; $arg_cnt++; }
}

# Hello!
print "\n\n\nYou've started the gene diversity script!\n\n";
my $dbname = "";

# Do we have everything we need to start?

# a table of loci and isolates and whether to transpose it or no
if(! defined $fTable) 
{ 
	print "Where is the table with isolates and loci?\n"; 
	$fTable = <STDIN>; 
	chomp $fTable; $fTable =~ s/\s+$//; # removes white spaces and line breaks
	
	
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
if(! -e $fTable)  { Usage("Input table file does not exist: $fTable"); exit; }


# a directory with FASTA files
if(! defined $dFAS)  
{ 
	print 	"\n\nWhere is the directory that includes the FASTA sequences?\n".
			"Alternatively, enter the name of the sequence database in BIGSDB ".
			"(for example: pubmlst_chlamydiales_seqdef) for exporting out of the BIGSdb API.\n";
	$dFAS = <STDIN>; 
	chomp $dFAS; $dFAS =~ s/\s+$//; # removes white spaces and line breaks
	
	if ( $dFAS =~ /\./ ) # if it looks like a file address
	{ if(! -e $dFAS)  { Usage("Input FASTA directory doesn't exist: $dFAS"); exit; } }
	else 
	{ $dbname = $dFAS; }
}

# a directory where results can go
if(! defined $dOut)  
{ 
	my @tableaddress = split ("/", $fTable);
	$tableaddress[-1] =~ s/\.[^.]+$//;
	$tableaddress[-1] = "GeneDiv-" . $tableaddress[-1];
	
	print 	"\n\nName the directory where you'd like the results to go.\n".
			"If you can't decide, leave this blank and a new folder, in the same place as your ".
			"table is held will be created, named '$tableaddress[-1]'\n";
	$dOut = <STDIN>;

	if ( $dOut =~ /\w/ )
	{ if(  -e $dOut)  { Usage("Output directory already exists: $dOut"); exit; } }
	else 
	{
		#push (@tableaddress, "/" ); 
		$dOut = join ("/", @tableaddress); 
	}

}

# Now that we've asked for everything, let's check that it's all what we expect... 
print "\n\n\nLet's check if everything is correct before we start...\n\n";

if ( $transpose eq "No" )
{ print	"Your table of isolates as columns and loci as rows is held at $fTable.\n"  }
elsif ( $transpose eq "Yes" )
{ print	"Your table of loci as columns and isolates as rows is held at $fTable.\n"  }


# including that folder we'll be putting all the results into
if( -e $dOut) { Usage("Output directory already exists: $dOut"); exit; }
elsif ( mkdir $dOut )
{ 
	print "Your results will be in the new directory $dOut\n"; 
}
else { Usage("Could not create output directory: $dOut"); exit; }

open ( RESULTS, '>', $dOut."/ResultsTable.txt" ) or die "$dOut /ResultsTable.txt"; # then also open where the reduced one will go
	print RESULTS join ("\t", ("Locus","Missing","Paralogous","Count-Nuc","Count-AA",) );

# and where the FASTA sequences are / will be
if ( $dbname =~ /\w/) # if given database name instead of folder, create that folder now
{ 
	print "FASTA sequences will be downloaded from the BIGSdb API from database $dbname.\n"; 
	
	if ( mkdir $dOut."/BIGSdb-FASTA/" )
	{ print "They'll be going into a folder called /BIGSdb-FASTA within your results folder. \n" }
	else 
	{ Usage("Could not create FASTA exports folder: $dOut/BIGSdb-FASTA/"); exit; }
}
else 
{ print "Your directory of FASTA sequences is at $dFAS.\n" }


# Give the option to get out if it goes in fact look dodgy
print "\nLeave blank to continue, enter any letters to exit.\n\n";
my $confirmation = <STDIN>;
if ( $confirmation =~ /\w/ )
{ Usage("You chose to exit: $confirmation"); exit;  }



# Actually starting the work now


# Reading the table in
my $tableref = ReadTableIn ( $fTable );
my @aoaTable = @$tableref; #puts table, in array of arrays format, into aoaTable

if ( $transpose eq "Yes" ) # transposes the table and puts it back in aoaTable
{
	$tableref = TransposeTable ( \@aoaTable );
	@aoaTable = @$tableref;
}


# Going through table

mkdir $dOut."/Observed-FASTA/" or die "Could not create filtered FASTA folder: $dOut/Observed-FASTA/";

foreach my $locusrow (@aoaTable) # loop per locus
{
	# empty hash to find unique values, values to be found
	my %unique_alleles = ();
	my $paralogous = 0; 
	my $missing = 0;
	my $countnuc = 0;

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

	# find the file where the original FASTA sequences are
	my $fullFAS; 
	
	if ( $dbname =~ /\w/ )
	{
		GetFASTASeqs ( $dbname, $locusname ) ;
		$fullFAS = $dOut."/BIGSdb-FASTA/".$locusname.".FAS";
	}
	else
	{	$fullFAS = $dFAS."/".$locusname.".FAS"; }
	
	# name where the file with the reduced number of FASTA sequences will go
	my $reducedFAS = $dOut."/Observed-FASTA/".$locusname.".FAS";


	if ( open(FULLFASTA, $fullFAS) ) # if can actually find the full FASTA file 
	{
		open(REDFASTA, '>', $reducedFAS) or die "Cannot open $reducedFAS\n"; # then also open where the reduced one will go
		{       
			my $save = 1; # whether to copy the sequence following a >identifier line
			
			while ( my $line = <FULLFASTA> ) # reading through original FASTA
			{
				if ( $line =~ /^>/ ) # if line is >identifier
				{
					chomp $line; 
					
					$line =~ /([^>_]+)$/; $line = $1; #removes > and anything before the last _
					
					if ( exists($unique_alleles{$line}) && $unique_alleles{$line} >= $dup ) # exists in hash of wanted loci
					{
						print REDFASTA "\n>", $line, "\n";
						$save = 1;
						$countnuc ++;
						$unique_alleles{$line} = 0; # set value to zero (normally would be undef or 1+) as check that has been read
					} 
					else # if not in the hash created from the table, indicating whether this allele is observed
					{ $save = 0; } 
				}
				
				elsif ( $line =~ /^([A-Za-z*-])/ ) # if is a sequence line (starts with a letter, - or *) 
				{
					chomp $line; 
					if ($save == 1) # copy only if "save" is turned on because last ID was wanted
					{ print REDFASTA $line; }
				}	
			}
						
			if ( exists($unique_alleles{"0"}) ) # move count away from empty if any missing alleles were seen
			{	$missing = $unique_alleles{"0"} } # gives the value in the frequency hash when the key is allele "0", the missing allele
			
		}  
		close(FULLFASTA);
		close(REDFASTA);
	
		# check if all alleles seen in table were copied 
		foreach my $key ( sort keys %unique_alleles ) 
		{ 
			if ( $unique_alleles{$key} >= $dup && $key ne 0 ) # if frequency is still above cut-off (copied ones are reset to 0) and isn't the missing allele, 0
			{ print "Did not find sequence for locus $locusname, allele $key.\n"; }
		}			
	}
	
	else # if there isn't a FASTA file to copy over, give warning plus example alleles (might be just "0" which would explain missing file)
	{ 
		print "$locusname did not exist as a FASTA file. Alleles in table were..."; 
		my $sampleallele = 0; 
		foreach my $key ( sort keys %unique_alleles ) 
		{ 
			if ( $sampleallele < 5 )
			{ print " $key"; $sampleallele++; }
			else # if already printed 5 allele numbers, just leave it
			{ print " etc."; last; }
		 }
		print "\n";
	}
	
	print RESULTS join ("\t", ($locusname, $missing, $paralogous, $countnuc) );
} # close per locus loop


# Translating 
mkdir $dOut."/Translated-FASTA/" or die "Cannot create /Translated-FASTA/ folder";

# opening the directory where the observed alleles are and getting all the file names
opendir (OBSDIR, $dOut."/Observed-FASTA/" ) or die "Cannot open directory: $!";
	my @files = readdir OBSDIR;
	@files = grep(/^([A-Za-z0-9])/,@files);
	if ($#files < 1)
	{ die "It looks like you have no FASTA files to align. If their file names don't begin with letters or numbers they are not being included.\n";}
closedir OBSDIR;


# Loop that translates all the files and puts them in new directory
print "Translated up to...\n";
my @lengths;

foreach my $file (@files)
{
 	my $infile  = $dOut . "/Observed-FASTA/"   . $file; # file with FASTA sequences in nucleotide form
 	my $outfile = $dOut . "/Translated-FASTA/" . $file; # file with FASTA sequences in amino acid form
 	
	my %unique = (); # empty hash to check for synonymous mutations 

	open(NUCLEOTIDE, $infile) 		or die "Cannot open $infile\n";
	open(AMINOACID, '>', $outfile) 	or die "Cannot open $outfile\n";

 	while ( my $line = <NUCLEOTIDE> )
	{
		if ( $line =~ /^>/ ) # if line is title of sequence
		{ 
			print AMINOACID $line;
		}
		elsif ( $line =~ /^[A-Z]/ ) # if line is sequence
		{
			chomp $line;
			my $peptide = Translate($line, "1");
			
			# for later checking the count of unique amino acids
			if ( !exists($unique{$peptide})) #if allele doesn't exist in unique-hash
			{ $unique{$peptide} = 1; } #then add it with allele number as key, value as 1
			
			print AMINOACID $peptide, "\n";
		}	
	}
	
	# counting number of unique amino acids and measuring the length of the locus
 	my ($locusname, $extension) = split (/\./, $file);
 	my $countaa = 0; 
 	my @lengths;
 	
 	foreach my $key ( sort keys %unique ) # goes through hash where keys are the amino acid sequences for each allele
 	{ 	
 		$countaa++;
 		push (@lengths, (3 * length($key) ) );
 	}

	if ( $countaa == 0 ) # if no alleles for that locus, have an escape route
	{	my ($min, $max, $avg) = 0, 0, 0; }
	else 
	{
		use List::Util qw( min max sum );
		my $min = min @lengths;
		my $max = max @lengths;
		my $avg = (sum @lengths) / $countaa;
 	}
 	
 	close NUCLEOTIDE; close AMINOACID;
 	
 	print RESULTS join ("\t", ($locusname, $missing, $paralogous, $countnuc) ); $count . "," . $avg . "\n";

 	#So you have something to watch while it runs... 
 	print "\r$locusname";
}

print "\nCompleted translation\n";





# closes per-locus loop

=begin GHOST 





# Getting the list of loci from the table
my $headerref = $aoaTable[0]; #dereferences the first line of the table
my @locuslist = (); 

foreach ( @$headerref ) # loops through each item in first row
{
	if ( /^(\S+)/ ) # takes only the first "word", and only if it exists in the %groups hash
	{ push @locuslist, $1; } # silently add name of group to array as in order
}




if ( $dbname =~ /\w/)
{
	print "We'll start by exporting those FASTA files from $dbname.\n";
	
	# first, make and go to the folder where they'll go
	
	if ( mkdir "BIGSdb-FASTA/" )
	{ print "They'll be going in a folder called BIGSdb-FASTA within your results folder. \n" }
	else 
	{ Usage("Could not create FASTA exports folder: $dOut/BIGSdb-FASTA/"); exit; }
	
	GetFASTASeqs ( $dbname, \@locuslist ) 
}


close(UNIQUENUC);

	print UNIQUENUC "locus,count-nuc,missing\n";

=cut GHOST

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
	my $locus	= $_[1];

	my $url = "http://rest.pubmlst.org/db/".$dbname."/loci/".$locus."/alleles_fasta";
	my $file = $dOut."/BIGSdb-FASTA/".$locus.".FAS";
		
	getstore($url, $file);
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

GeneDiv

Description:
	  
Usage:
gene-diversity.pl [ options ]

Examines the diversity of a set of genomes across tagged loci. 

-tin: Table as tab-separated, with loci and isolates as rows and columns, or reverse if "transpose" option used.
-transpose: Term included if the table has loci as columns and isolates as rows.  
-FASTA: Directory with FASTA files where "locusname.FAS" is the file name, like BACT000001.FAS
-dout: Directory where all results will be saved.
-dup: Value for the mininum frequency for allele appearance to be included. Default is "1".
EOU

print "\n\nQuit because: $message\n\n";
print "ARGV was " . join (", ", @ARGV) . "\n";
}
