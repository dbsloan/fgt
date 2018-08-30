# fgt
Perl script to implement the four-gamete test

The accompanying FGT.pl script reads in a sequence alignment file in Fasta format, identifies all parsimony-informative sites, and performs a standard four-gamete test between all pairwise combinations of identified sites by counting the number of genotypic combinations found in the dataset for each site pair. Individual sites with more than two alleles are ignored.

## Requirements

The FGT.pl script is written in Perl and requires the sloan.pm module (https://github.com/dbsloan/perl_modules).

## Usage

The script can be called from a command-line terminal session as follows:

perl FGT.pl input_fasta_file > output_file 2> log_file

The output_file is a tab-delimited text file summarizing the results of the four-gamete test for each site pair. The log_file reports some basic statistics about the number of parsimony-informative sites that were found, as well as any warnings or errors.