#!/usr/bin/perl
use strict;

print 	"Please give me a tab-delimited file to transpose.\n";
my $infile = <STDIN>; #asks for infile
chomp ($infile);
my $outfile = $infile . "-t"; #names outfile 

print "You've put in $infile. You'll get results at $outfile.\n";

if (-e $outfile) #if outfile already exists
{ die("The results file $outfile already exists. Delete, rename or move it.\n"); }

my @original = (); #where data will go in the beginning
my @transposed = (); #where data will go after being transposed
my $columncount = 0;
my $rowcount = 0;

open(INFILE, $infile) or die "Cannot open $infile\n";
	while ( my $line = <INFILE> )		
	{
        chomp($line);
        $line =~ s/\r\n?|\n//g; #just in case, removes all other types of line break
        $line =~ s/,/\t/g; #if in comma format, switches it to tabs
        
        my @row = split ('\t', $line);
        $rowcount = @row;
        $columncount ++;
        push ( @original, [ @row ] );
	} 
close(INFILE);

for my $row (@original) #actually does the transposition
{
 	for my $column (0 .. $rowcount) 
 	{
   		push(@{$transposed[$column]}, $row->[$column]);
  	}
}

@transposed = splice (@transposed, 0, $rowcount); #removes empty rows if more isolates than loci

open(OUTFILE, '>>', $outfile) or die "Cannot open $outfile\n"; #opens your results file as appending to it
	for my $new_row (@transposed) 
	{
		for my $new_col (@{$new_row}) 
		{
			print OUTFILE $new_col, "\t";
		}
		print OUTFILE "\n";
	}
close(OUTFILE);

#largely borrowed from http://stackoverflow.com/questions/3249508/transpose-in-perl