## Options which may be added to command line:
* "-help": Will exit the script and call up the Usage instructions  
* "-table": A file containing a tab-separate table of loci vs. isolates with allele IDs in the cells  
* "-transpose": "Yes" if the table has loci as columns and isolates as rows, "No" otherwise.  
* "-FASTA": Directory with complete FASTA files where "locusname.FAS" is the file name, like BACT000001.FAS  
* "-dbname": If grabbing FASTA sequences from BIGSdb, the name of the database, usu. pubmlst\_*genusname*\_seqdef  
* "-dboption")	{ $FASTAoption 	= $ARGV[$i+1] || ''; $arg_cnt++; }
* "-out")			{ $dOut  		= $ARGV[$i+1] || ''; $arg_cnt++; }
* "-dup")			{ $dup  		= $ARGV[$i+1] || ''; $arg_cnt++; }
* "-mafft")		{ push (@mafftarg, $ARGV[$i+1]) ; $arg_cnt++; }


## Parameters in the Results table:
* *Locus*: name of locus as seen on input table  
* *Missing*: number of isolates for which no allele ID was present (missing, untagged, incomplete, etc.)  
* *Paralogous*: number of isolates in which the locus had at least two tagged alleles  
* *CountNuc*: count of unique nucleotide sequences for specified locus  
* *CountAA*: count of unique amino acid sequences for specified locus  
* *MinLength*: length in nucleotide letters of shortest nucleotide sequence  
* *MaxLength*: length in nucleotide letters of longest nucleotide sequence  
* *AvgLength*: average length in nucleotide letters of unique nucleotide sequences (used in Calculations)  
* *NonVarNuc*: number of sites in nucleotide alignment were no variation was seen (marked "\*" in ClustalW format)  
* *NonVarAA*: number of sites in amino acid alignment were no variation was seen (marked "\*" in ClustalW format)  


## Parameters added in the Calculations table:
* *AllelicDiv*: number of unique nucleotide sequences found divided by average length (greater if more diverse)  
* *ADivNM*: the previous parameter subtracted by the genome-wide average 
* *VSitesNuc*: proportion of sites in aligned nucleotide sequences which had some variation  
* *VSitesAA*: proportion of sites in aligned nucleotide sequences which had some variation  
* *RatioCount*: ratio of unique nucleotide to unique amino acid sequences (increases with diversifying selective pressure)  
* *RatioVS*: ratio of proportions of unique sites in nucleotide to amino acid sequences (increases with diversifying selective pressure)  


## Possible additions:
* measure of GC content  
* accepting loci categories  
* option of checking for variable sites by \*, :, or .  


## Example
