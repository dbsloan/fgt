#!/usr/bin/perl

use strict;
use warnings;
use List::Util qw (min);
use sloan;

my $usage = "\nUSAGE: perl $0 input_fasta > output 2> log_file\n\n";

my $file = shift or die ($usage);

my %fasta = fasta2hash($file);

my $length = 0;
foreach (keys %fasta){
	if ($length){
		length ($fasta{$_}) == $length or die ("\nERROR: fasta input sequences do not have identical lengths.\n\n");
	}else{
		$length = length($fasta{$_});
	}
}

my %sitesHoH;

my %multinuc;
my %singletons;

for (my $i = 0; $i < $length; ++$i){

	my %nuc_hash;
	
	foreach (keys %fasta){
		++$nuc_hash{uc(substr($fasta{$_}, $i, 1))};
	}

	scalar (keys %nuc_hash) == 1 and next;
	scalar (keys %nuc_hash) > 2 and $multinuc{$i+1} = 1 and next;
	if (scalar (keys %nuc_hash) == 2){
		if (min (values %nuc_hash) == 1){
			$singletons{$i+1} = 1;
		}else{
			foreach my $header (keys %fasta){
				$sitesHoH{$i+1}->{$header} = uc(substr($fasta{$header}, $i, 1));
			}
		}
	}else{
		die ("\nERROR: nucleotide counting failed at site $i (0 indexed)\n\n");
	}
}

my $singleton_count = scalar (keys %singletons);
my $multinuc_count = scalar (keys %multinuc);
my $parsimony_count = scalar (keys %sitesHoH);
my $polymorphic_count =  $singleton_count + $multinuc_count + $parsimony_count;


print STDERR "Input File: $file\n\nAlignment Length: $length\nPolymorphic Sites: $polymorphic_count\nSingleton sites: $singleton_count\nExcluded sites with >2 alleles: $multinuc_count\n\nSites analyzed for FGT: $parsimony_count\n";


print "Site1\tSite2";

foreach (sort keys %fasta){
	print "\t$_";
}

print "\tAlleleCount\n";

my @pos_array = sort {$a <=> $b} keys %sitesHoH;

for (my $i = 0; $i < scalar (@pos_array); ++$i){
	for (my $j = $i+1; $j < scalar (@pos_array); ++$j){
		print "$pos_array[$i]\t$pos_array[$j]";
		my %dinuc_hash;
		foreach (sort keys %fasta){
			my $dinuc = uc($sitesHoH{$pos_array[$i]}->{$_}) . uc($sitesHoH{$pos_array[$j]}->{$_});
			++$dinuc_hash{$dinuc};
			print "\t$dinuc";
		}
		my $allele_count = scalar (keys %dinuc_hash);
		print "\t$allele_count\n";
	}
}


