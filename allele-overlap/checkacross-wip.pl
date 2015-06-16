#!/usr/bin/perl
use strict;

#asks for files, removes linebreak and white spaces, puts each file as one item in an @input
my $statement = <<END; 
Please list, separated by commas, the files you want to compare.
Each file should contain one group's list of unique alleles, in "BACT000001,12,5345,7424" format.
END
print $statement; 
my $infiles = <STDIN>; 
chomp $infiles;
$infiles =~ s/\s//g; 
my @input = split(/,/, $infiles); 

#asks for results name, removes linebreak, checks if already exists, if yes error
print "What do you want to call the results?\n";
my $outfile = <STDIN>;
chomp $outfile;
if (-e $outfile) 
{ die("The results file $outfile already exists. Delete, rename or move it.\n"); }

#checks if at least two files were entered, tells you what and how many, where results will be
if (@input < 2)
{ die ("You need at least two files to compare."); }
my $filecount = $#input + 1;
print "You are comparing ".join (", ", @input).", $filecount in total.\nYour results will be at $outfile.\n";

#creates empty hash table for file information to go, and two empty results arrays
our %locitable = ();
my @results = (); 
my @sharedresults = (); 

#goes through each file, puts their information into %locitable
# %loci table keys are locus names, values are allele numbers (all groups' alleles together)
foreach my $file (@input)
{ populatesHashfromInfiles ($file, \%locitable); }

print %locitable; 

#goes through each array (allele numbers) in hash of arrays (list of loci)
foreach my $locusName (keys %locitable)
{
	my %allele_frequencies = (); #creates empty unique-hash
	my @sharedalleles = (); #creates empty shared allele array
	
	my @alleleArray = @{ $locitable{$locusName} }; #puts each allele in that locus as one item in an array
	
	foreach my $allele (@alleleArray) #go through each allele in locus
    {
    	if ( !exists($allele_frequencies{$allele})) #if allele doesn't exist in unique-hash
    	{ $allele_frequencies{$allele} = 1; } #then add it with allele number as key, value as 1
    	else
    	{ 
    		$allele_frequencies{$allele} ++; #otherwise up the count of the allele
    		if ($allele_frequencies{$allele} == 2) { 
    			push @sharedalleles, $allele; 
    		} #when allele is known to be shared, add to @sharedallele
    	} 
	}
	
	#takes out all the allele frequencies (values) from the allele_frequencies hash, sorts them into @frequencies
	my @frequencies = sort {$a <=> $b} values(%allele_frequencies);
	
	my $totalAlleleCount = @frequencies; #count of items in @frequencies is the count of unique alleles across all groups
	my $specificAlleleCount = $totalAlleleCount - @sharedalleles; #count of specific alleles is all minus shared
	my $sharedAlleleCount = @sharedalleles; #count of shared alleles
	
	if ($#frequencies == 0) #if the count of allele frequencies is one
	{ push @results, "$locusName: no variation\n"; } #then that mean there is only one allele in all groups
	
	else
	{
		my $first = shift(@frequencies); #first value in frequencies into $first
		my $last = pop(@frequencies); #last value in frequencies into $last
		
		if ($first == $last) #if first and last are all the same, then all are the same since array is sorted
		{
			if ($first == 1) #if all alleles occur only once across all groups
			{ push @results, "$locusName: all alleles group specific, 0\n"; } #then must be specific to the group
			elsif ($first == $filecount) #if all alleles occur in each group once
			{ push @results, "$locusName: all alleles shared among all groups, 1\n"; } #then must all overlap
		}
		else #mix of frequencies
		{
			if ($last == $filecount) #if highest frequency is total number of files compared
			{ push @results, "$locusName: some alleles shared among all groups, $sharedAlleleCount/$totalAlleleCount\n"; }
			else #no alleles shared among all of the groups
			{ push @results, "$locusName: some alleles shared among some groups, $sharedAlleleCount/$totalAlleleCount \n"; }
			push @sharedresults, "$locusName: shared alleles include @sharedalleles\n";
		}
	}
}

@results = sort {$a cmp $b} @results; #sorts results by locus name
@sharedresults = sort {$a cmp $b} @sharedresults; #sorts results by locus name

open(OUTFILE, '>>', $outfile); #open results file
	print OUTFILE "Results for " . join (', ', @input) . "\n\n";
	print OUTFILE @results; #prints results array to file
	print OUTFILE "\n\n\n";
	print OUTFILE @sharedresults;
close(OUTFILE); 


sub populatesHashfromInfiles
{
	my $infile = @_[0];
	my $hash_ref = @_[1];
	
	print "Infile is $infile, hash is %hash.\n";
	
	open(INFILE, $infile) or die "Cannot open $infile\n";
		while ( my $line = <INFILE> )	{        
			#cleans up line breaks, makes sure in comma format
			$line =~ s/\r\n?|\n//g;
      		$line =~ s/\t/,/g;
      		
      		#separates this line of the file into multiple values in @linedata array
			my @linedata = split(/,/, $line); 

			#removes the header and places it into $header for safekeeping
			my $header = shift(@linedata);  
			
			#looks through the rest of the array, putting each value in %hash with header as key
			foreach my $value (@linedata) {
				push @{%$hash_ref{$header}}, $value;
			}	
		}  
	close(INFILE);
}