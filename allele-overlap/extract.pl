#!/usr/bin/perl
use strict;

print 	"Please give me a csv file with isolates on columns, loci on rows, and allele numbers as the values.\n";
my $infile = <STDIN>; #asks for infile
chomp ($infile);
my $outfile = $infile . "-unique"; #names outfile

print "You've put in $infile. You'll get results at $outfile.\n";

if (-e $outfile) #if outfile already exists
{ die("The results file $outfile already exists. Delete, rename or move it.\n"); }

open(INFILE, $infile) or die "Cannot open $infile\n";
	<INFILE>; #ignores header of file (usually "Locus", "45045 | ERR09342093", etc.)
	while ( my $line = <INFILE> ) #goes through each row in the file at a time			
	{
	    chomp($line); #removes line break from end of line
		$line =~ s/\r\n?|\n//g; #just in case, removes all other types of line break
        $line =~ s/\t/,/g; #if in tab format, switches it to commas
        
        my @linedata = grep {length()} split(/,/, $line); #put each allele into one item of an array, skips empty cells 

        my $locus = shift(@linedata); #takes out first array item (BACTXXXX) puts it into $locus scalar
		
    	my %unique_alleles = (); #empties unique-counting hash table
    	
    	foreach my $allele (@linedata) #goes through each allele in array
    	{
    		if ( !exists($unique_alleles{$allele})) #if it's not in the hash table
    		{ $unique_alleles{$allele} = 1; } #creates it! 
		}
		
		my @unique = (); #empties unique array
		
		foreach my $key ( sort keys %unique_alleles ) #goes through each item in hash table
		{
  			push (@unique, $key) #puts hash keys into array, so puts each unique allele number into array
		}
		
		unshift (@unique, $locus); #puts BACTXXXX back in front of array so you know which locus the alleles are at
		
		open(OUTFILE, '>>', $outfile) or die "Cannot open $outfile\n"; #opens your results file as appending to it
			print OUTFILE join (',',@unique) . "\n"; #print results to outfile, as "BACTXXXX,12,4234,234"
		close(OUTFILE);
		
	} 
close(INFILE);

print 	"All done!\nCheck your results at $outfile. Each line begins with the locus name and contains all unique alleles at that locus.\n";
