#!/usr/bin/perl -w

use strict;
use warnings;
use Carp;
use lib "lib";
use Weitz::BioDiversity::Util;
use Weitz::BioDiversity::Time;

my $param_file = "params.conf";

my %params = Weitz::BioDiversity::Util->readParamFile($param_file);

if($params{"download"} eq "Y" || $params{"download"} eq "y"){
	print "\n*****************************************************************\n";
	print "Downloading data. It will take a while. Thanks for your patience.\n";	
	Weitz::BioDiversity::Util->downloadViralGenomes();
	print "Download complete.\n";
	print "*****************************************************************\n";		
}

my @genome_file_names = Weitz::BioDiversity::Util->readGenomesFileNames();

my %all_genome_length = Weitz::BioDiversity::Util->getAllViralGenomeLength(\@genome_file_names);

# Get the keys of the above hash sorted in descending order of the genome length
my @sorted_genomes = Weitz::BioDiversity::Util->getSortedKeysOfHashByDescValue(%all_genome_length);
if($params{"class"} > @genome_file_names){
	croak "\n\nCommunity size can not exceed total number of genomes.\nSelect a value <= ".@genome_file_names.".\n\n";
}

my $comm_aref = Weitz::BioDiversity::Util->generateCommunity(\%params);
my %sample = Weitz::BioDiversity::Util->generateSample($params{"sample"}, $comm_aref);

my $end_idx = $params{"class"} - 1;
my @selected_genomes = @sorted_genomes[0..$end_idx];

my $sum = 0;
my $sample_aref = $sample{"sample"};
$sum += $_ for @$sample_aref;
print "\n*****************************************************************\n";
print "Simulating ". $sum." reads from ".@selected_genomes." genomes, please be patient.\n";
my $start_time = Weitz::BioDiversity::Time->new()->getTime();
print "Start time: ", Weitz::BioDiversity::Time->new()->getLocalTime(), " \n";
my %sim_results = Weitz::BioDiversity::Util->simulateReads(\@selected_genomes, $sample{"sample"}, $params{"read"});
print "End Time: ", Weitz::BioDiversity::Time->new()->getLocalTime(), " \n";	
my $end_time = Weitz::BioDiversity::Time->new()->getTime();
my $run_time = $end_time - $start_time;
print "Read generation run time:: $run_time milliseconds \n\n";	
print $sim_results{"read_count"}." reads generated.\n";
print "*****************************************************************\n";

my $blast_db = Weitz::BioDiversity::Util->makeBlastDB($sim_results{"read_file"});
print "\n*****************************************************************\n";
print "Running pairwise blast for ", $sim_results{"read_count"}, " reads. It will take a while.\n";
print "Please be patient.\n";
my $s_time = Weitz::BioDiversity::Time->new()->getTime();
print "Blast started at ", Weitz::BioDiversity::Time->new()->getLocalTime(), " \n";
my $blast_out_file = Weitz::BioDiversity::Util->runBlast($blast_db, $sim_results{"read_file"}, $sim_results{"read_count"}, $params{"blast"});
print "Blast ended at ", Weitz::BioDiversity::Time->new()->getLocalTime(), " \n";
my $e_time = Weitz::BioDiversity::Time->new()->getTime();
my $r_time = $e_time - $s_time;
print "BLAST run time:: $r_time milliseconds \n\n";
print "*****************************************************************\n";

my %blast_result = Weitz::BioDiversity::Util->getBlastResult($blast_out_file);

print "Distribution of reads in sample:\n@$sample_aref\n\n";
print $blast_result{"TM"}." unique pairwise matches found.\n";
print $blast_result{"FP"}." false positive matches found.\n\n";
my %false_values = Weitz::BioDiversity::Util->estimateFalseValues(\%blast_result, $sample{"sample"}, $sim_results{"read_count"});
print "For Evalue of ".$params{"blast"}->{"evalue"}.", false +ve rate is ".$false_values{"positive"}." and false -ve rate is ".$false_values{"negative"}."\n";
my %similarity = Weitz::BioDiversity::Util->estimateSimilarity($sample{"community"}, $blast_result{"TM"}, $sim_results{"read_count"});
my $correct_ss;
if($false_values{"negative"} < 1){
	$correct_ss = sprintf("%.6f",$similarity{"SS"}*(1/(1-$false_values{"negative"})));
}else{
	$correct_ss = "1/0 ";
}
print "\n*****************************************************************\n";
print "Community Similarity			:: ".$similarity{"CS"}."\n";
print "Sample Similarity   			:: ".$similarity{"SS"}."\n";
print "Corected Sample Similarity		:: ".$correct_ss."\n";
print "\n*****************************************************************\n";

#Weitz::BioDiversity::Util->cleanRunDir($blast_out_file);
