#!/usr/bin/perl -w

use strict;
use warnings;
use Carp;
use lib "lib";
use Weitz::BioDiversity::Util;
use Weitz::BioDiversity::Time;
use Getopt::Long;
use File::Basename;
use POSIX;

# receive command line arguments to run the simulation
my($max_iterations, $min_evalue, $max_evalue, $min_reads, $max_reads, $read_increment, $start_index, $end_index, $help) = (2,0.000001, 0.001, 200, 1000, 200, 0, 2, 0);

&processCMDLine;

my $param_file = "params.conf";

my %params = Weitz::BioDiversity::Util->readParamFile($param_file);

if($params{"download"} eq "Y" || $params{"download"} eq "y"){
	Weitz::BioDiversity::Util->downloadViralGenomes();
}

my @genome_file_names = Weitz::BioDiversity::Util->readGenomesFileNames();

my %all_genome_length = Weitz::BioDiversity::Util->getAllViralGenomeLength(\@genome_file_names);

# Get the keys of the above hash sorted in descending order of the genome length
my @sorted_genomes = Weitz::BioDiversity::Util->getSortedKeysOfHashByDescValue(%all_genome_length);

if($params{"class"} > @genome_file_names){
	croak "\n\nCommunity size can not exceed total number of genomes.\nSelect a value <= ".@genome_file_names.".\n\n";
}
if((($end_index - $start_index) + 1) > @genome_file_names){
	croak "\n\nCommunity size can not exceed total number of genomes.\nSelect a value <= ".@genome_file_names.".\n\n";	
}

# Simulate over multiple read sizes , evalues 
#&simOneDisManyReadsAndEvalues;

# Simulate over different community sizes over a fixed evalue and number of reads
&simManyCommunity;

# Simulate with a fixed distribution
#&simFixedDistribution;

sub simFixedDistribution{
	my $c_size = (($end_index - $start_index) + 1);
	if($c_size != 2){
		die "Invalid start and end index. The difference should be equal to 2\n";
	}
	my $output = "Distribution Fixed: \n";
	my $size = $params{"sample"}; 
	print "\n\nRunning Simulation for a distribution with fixed sample size of $size\n\n";
	print "Start time: ", Weitz::BioDiversity::Time->new()->getLocalTime(), " \n\n";

	for(my $g = 1; $g <= 9; $g++){
		# genomes should be sliced from selected from the 
		my @selected_genomes = @sorted_genomes[$start_index..$end_index];
		$params{"class"} = (($end_index - $start_index) + 1);
		my $comm_aref = getFixedDistribution($g, 10-$g);

		$output = $output."Distribution probabilities: @$comm_aref \n\n";
		my %sample;
		print "\n\nRunning $max_iterations simulations for evalue of ".$params{"blast"}->{"evalue"}." on ". @selected_genomes." communities with ". $size. " reads \n\n";
		%sample = Weitz::BioDiversity::Util->generateSample($size, $comm_aref);
		#prepare a string to be written to a file.
		$output .= buildHeader($size, \@selected_genomes, $sample{"sample"});
		$output = $output.$params{"blast"}->{"evalue"};
		print "Evalue ".$params{"blast"}->{"evalue"}." :: Iteration: ";
		for(my $i = 0; $i < $max_iterations; $i++){
			print "$i,";
			my %sim_results = Weitz::BioDiversity::Util->simulateReads(\@selected_genomes, $sample{"sample"}, $params{"read"});
			my $blast_db = Weitz::BioDiversity::Util->makeBlastDB($sim_results{"read_file"});
			my $blast_out_file = Weitz::BioDiversity::Util->runBlast($blast_db, $sim_results{"read_file"}, $sim_results{"read_count"}, $params{"blast"});
			my %blast_result = Weitz::BioDiversity::Util->getBlastResult($blast_out_file);
			my %false_values = Weitz::BioDiversity::Util->estimateFalseValues(\%blast_result, $sample{"sample"}, $sim_results{"read_count"});
			my %similarity = Weitz::BioDiversity::Util->estimateSimilarity($sample{"community"}, $blast_result{"TM"}, $sim_results{"read_count"});
			my $correct_ss;
			if($false_values{"negative"} < 1){
				$correct_ss = sprintf("%.6f",$similarity{"SS"}*(1/(1-$false_values{"negative"})));
			}else{
				$correct_ss = "1/0 ";
			}
			$output = $output."\t".$false_values{"positive"}.",".$false_values{"negative"}.",".$similarity{"CS"}.",".$similarity{"SS"}.",".$correct_ss;
			Weitz::BioDiversity::Util->cleanRunDir($blast_out_file);
		}
		print "\n";
		$output = $output."\n";
		$output .= "\n\n*******************************************************************************\n\n";
	}
	&writeResults($output);
}


sub simManyCommunity{
	my $output = "Distribution : ".$params{"distribution"}."\n";
	my $size = $params{"sample"}; 
	print "\n\nRunning Simulation for a distribution with fixed sample size of $size\n\n";
	print "Start time: ", Weitz::BioDiversity::Time->new()->getLocalTime(), " \n\n";

	for(my $g = 10; $g < @sorted_genomes;){
		# genomes should be sliced from the sorted list 
		my @selected_genomes = @sorted_genomes[0..$g-1];
		$params{"class"} = $g;
		my $comm_aref = Weitz::BioDiversity::Util->generateCommunity(\%params);

		$output = $output."Distribution probabilities: @$comm_aref \n\n";
		my %sample;
		print "\n\nRunning $max_iterations simulations for evalue of ".$params{"blast"}->{"evalue"}." on ". @selected_genomes." communities with ". $size. " reads \n\n";
		%sample = Weitz::BioDiversity::Util->generateSample($size, $comm_aref);
		#prepare a string to be written to a file.
		$output .= buildHeader($size, \@selected_genomes, $sample{"sample"});
		$output = $output.$params{"blast"}->{"evalue"};
		print "Evalue ".$params{"blast"}->{"evalue"}." :: Iteration: ";
		for(my $i = 0; $i < $max_iterations; $i++){
			print "$i,";
			my %sim_results = Weitz::BioDiversity::Util->simulateReads(\@selected_genomes, $sample{"sample"}, $params{"read"});
			my $blast_db = Weitz::BioDiversity::Util->makeBlastDB($sim_results{"read_file"});
			my $blast_out_file = Weitz::BioDiversity::Util->runBlast($blast_db, $sim_results{"read_file"}, $sim_results{"read_count"}, $params{"blast"});
			my %blast_result = Weitz::BioDiversity::Util->getBlastResult($blast_out_file);
			my %false_values = Weitz::BioDiversity::Util->estimateFalseValues(\%blast_result, $sample{"sample"}, $sim_results{"read_count"});
			my %similarity = Weitz::BioDiversity::Util->estimateSimilarity($sample{"community"}, $blast_result{"TM"}, $sim_results{"read_count"});
			my $correct_ss;
			if($false_values{"negative"} < 1){
				$correct_ss = sprintf("%.6f",$similarity{"SS"}*(1/(1-$false_values{"negative"})));
			}else{
				$correct_ss = "1/0 ";
			}
			$output = $output."\t".$false_values{"positive"}.",".$false_values{"negative"}.",".$similarity{"CS"}.",".$similarity{"SS"}.",".$correct_ss;
			Weitz::BioDiversity::Util->cleanRunDir($blast_out_file);
		}
		print "\n";
		$output = $output."\n";
		$output .= "\n\n*******************************************************************************\n\n";
		$g += 100;
	}
	&writeResults($output);
}


sub simOneDisManyReadsAndEvalues{
	# genomes should be sliced from selected from the 
	my @selected_genomes = @sorted_genomes[$start_index..$end_index];
	
	#Update the param{'class'} based on the value received on the commandline
	$params{"class"} = (($end_index - $start_index) + 1);
	
	my $comm_aref = Weitz::BioDiversity::Util->generateCommunity(\%params);
	my $output = "Distribution : ".$params{"distribution"}."\n";
	$output = $output."Distribution probabilities: @$comm_aref \n\n";
	my %sample;
	
	print "\n\nRunning Simulation for a distribution with different sample sizes\n\n";
	print "Start time: ", Weitz::BioDiversity::Time->new()->getLocalTime(), " \n\n";
	
	#my $size = $params{"sample"}; This is value is overridden by command line argument
	
	for(my $j = $min_reads; $j <= $max_reads;){
		print "\n\nRunning $max_iterations simulation for each evalue of ($min_evalue to $max_evalue) on ". @selected_genomes." communities with ". $j. " reads \n\n";
		%sample = Weitz::BioDiversity::Util->generateSample($j, $comm_aref);
		#prepare a string to be written to a file.
		$output .= buildHeader($j, \@selected_genomes, $sample{"sample"});
		
		for(my $e = $min_evalue; $e <= $max_evalue;){
			$output = $output.$e;
			print "Evalue $e :: Iteration: ";
			for(my $i = 0; $i < $max_iterations; $i++){
				print "$i,";
				my $blast_href = $params{"blast"};
				$blast_href->{"evalue"} = $e;
				my %sim_results = Weitz::BioDiversity::Util->simulateReads(\@selected_genomes, $sample{"sample"}, $params{"read"});
		
				my $blast_db = Weitz::BioDiversity::Util->makeBlastDB($sim_results{"read_file"});
		
				my $blast_out_file = Weitz::BioDiversity::Util->runBlast($blast_db, $sim_results{"read_file"}, $sim_results{"read_count"}, $params{"blast"});
				my %blast_result = Weitz::BioDiversity::Util->getBlastResult($blast_out_file);
				my %false_values = Weitz::BioDiversity::Util->estimateFalseValues(\%blast_result, $sample{"sample"}, $sim_results{"read_count"});
				my %similarity = Weitz::BioDiversity::Util->estimateSimilarity($sample{"community"}, $blast_result{"TM"}, $sim_results{"read_count"});
				my $correct_ss;
				if($false_values{"negative"} < 1){
					$correct_ss = sprintf("%.6f",$similarity{"SS"}*(1/(1-$false_values{"negative"})));
				}else{
					$correct_ss = "1/0 ";
				}
				$output = $output."\t".$false_values{"positive"}.",".$false_values{"negative"}.",".$similarity{"CS"}.",".$similarity{"SS"}.",".$correct_ss;
				Weitz::BioDiversity::Util->cleanRunDir($blast_out_file);
			}
			print "\n";
			$output = $output."\n";
			if($e < 1.0){
				$e *= 10;
			}else{
				$e += 1;
			}
		}
		$j += $read_increment;
		$output .= "\n\n*******************************************************************************\n\n";
	}	
	&writeResults($output);
}


sub writeResults{
	my($output) = @_;
	# write the genome length into a file
	my $outFile = "results/simulation_results_".Weitz::BioDiversity::Time->new()->toString().".txt";
	unless(open(OutFileHandle, ">$outFile")){
		croak("Failed to create $outFile file \n");
	}
	unless(print OutFileHandle $output){
		croak("Failed to write the content to $outFile file \n");
	}
	print "\n\nEnd time: ", Weitz::BioDiversity::Time->new()->getLocalTime(), "\n";
	print "Simulation results are written to the file $outFile \n";
}


sub buildHeader{
	my ($size,$gen_aref,$sample_aref) = @_;
	my @genomes = @$gen_aref;
	
	my $output = "";
	$output = $output."Total Sample Size: ".$size."\n\n";
	$output = $output."Genomes      : @genomes\n\n";
	$output = $output."Genome Length: ";
	foreach my $genome (@genomes){
		$output = $output."\t".$all_genome_length{$genome};
	}
	$output = $output."\n\n";
	$output = $output."Sample Size  : ";
	foreach my $read (@$sample_aref){
		$output = $output."\t".$read;
	}
	$output = $output."\n\n";
	$output = $output."Average Read Length: ".$params{"read"}."\n\n";
	$output = $output."$max_iterations iterations over a given Evalue\n\n";
	$output = $output."For each Iteration there are 5 values i.e. False Positive, False Negative, Community Similarity, Sample Similarity, Corrected Sample Similarity\n\n";	
	
	return $output;
}

sub processCMDLine{
	my $pFile = fileparse($0);
	if(@ARGV < 1){
		print "Invalid number of arguments. Check the usage: \n";
		&help($pFile); 
		exit;
	}
	
	GetOptions(
		'itr=i' => \$max_iterations,
	  	'm_eval=f' => \$min_evalue,
 	  	'x_eval=f' => \$max_evalue,
 	  	'm_read=i' => \$min_reads,
 	  	'x_read=i' => \$max_reads,
 	  	'r_incr=i' => \$read_increment,
 	  	's_idx=i' => \$start_index,
 	  	'e_idx=i' => \$end_index,
 	  	'help!' => \$help,
	) or die &help($pFile);
	
	#change the man and max number of reads based on the number of genomes selected. 
#	$min_reads = (($end_index - $start_index) + 1) * 50;
#	$max_reads = (($end_index - $start_index) + 1) * 200;
#	$read_increment = ceil($max_reads / 5);
	
	if($help){ &help($pFile); }
	
}

sub help{
	my($cmd_file) = @_;
	print "\n\nUSAGE:\n\n";
	print $cmd_file." [options]\n\n";
	print "Options are: \n\n";
	print "[-itr		value]			Total number of simulations to run\n";
	print "[-m_eval	value]			Minimum evalue to consider for each simulation\n";
	print "[-x_eval	value]			Maximum evalue to consider for each simulation\n";
	print "[-m_read	value]			Minimum size of reads to be ditributed across genomes\n";
	print "[-x_read	value]			Maximum size of reads to be ditributed across genomes\n";
	print "[-r_incr	value]			Number of reads to be incremented in each iteration\n";
	print "[-s_idx		value]			Start index of the genome list\n";
	print "[-e_idx		value]			End index of the genome list\n";
	print "[-help] [-h]				List all options available\n\n";
	print "Example:: $cmd_file -itr 10\n\n";
	exit;	
}

sub getFixedDistribution{
	my($a, $b) = @_;
	my $c = $a + $b;
	my @community = ($a/$c, $b/$c);
	return \@community;
}