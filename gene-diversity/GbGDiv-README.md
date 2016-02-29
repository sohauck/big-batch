## What is necessary for these scripts to work correctly
* Perl, including List::Util and LWP::Simple, download from https://www.perl.org/get.html
* MAFFT — tested on v7.221, download from http://mafft.cbrc.jp/alignment/software/ 
* R — tested on 3.0.2, download from https://cran.r-project.org/
  * ggplot2 — run 'install.packages("ggplot2")' from within R


## Options which may be added to command line:
* *-help*: Will exit the script and call up the Usage instructions  
* *-table*: A file containing a tab-separate table of loci vs. isolates with allele IDs in the cells  
* *-transpose*: 'Yes' if the table has loci as columns and isolates as rows, 'No' otherwise.  
* *-FASTA*: Directory with complete FASTA files where *locusname.FAS* is the file name, like BACT000001.FAS  
* *-dbname*: If grabbing FASTA sequences from BIGSdb, the name of the database, usu. pubmlst\_*genusname*\_seqdef  
* *-dboption*: '1' if complete FASTA directory exists, '2' if creating the same, '3' if pulling just what is needed  
* *-out*: Directory where all results will be saved.
* *-dup*: Value for the mininum frequency for allele appearance to be included. Default is '1'.
* *-mafft*: Can be used to add more parameters to MAFFT. 


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
* variable exclusion parameters for R step (currently cutoff at 10% of isolates tagged)
 

## Example
From the folder where the Perl '' and R '' scripts are, along with the 'examples' folder, run this: 
```
perl GbGDiv.pl -table examples/MYCO-table.txt -transpose Y -out examples/Trial-1/ -dbname pubmlst_mycobacteria_seqdef -dboption 2 -mafft --auto 
```
The results should be identical to Trial-0 folder.
