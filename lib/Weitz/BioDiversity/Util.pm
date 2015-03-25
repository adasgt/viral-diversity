package Weitz::BioDiversity::Util;

use strict;
use warnings;
use Carp;
use Weitz::BioDiversity::Powerlaw;
use Weitz::BioDiversity::Lognormal;
use Weitz::BioDiversity::Exponential;
use Weitz::BioDiversity::SampleCommunity;
use Weitz::BioDiversity::DistributionFactory;
use Weitz::BioDiversity::MetaSim;
use Weitz::BioDiversity::BlastPlus;
use Weitz::BioDiversity::Time;
use Weitz::BioDiversity::SimpleReadGenerator;

# Input to readParamFile is a parameter file.
# This method reads a parameter file and returns a key:value pair hash
sub readParamFile{
	my($class, $param_file) = @_;
	my %parameters = ();
	my %blast_params = ();
	unless(open(FileHandle, $param_file)){
		croak("Cannot open file $param_file. Make sure there is a file named $param_file.");
	}
	my @data = <FileHandle>;
	close(FileHandle);
	foreach my $line (@data){
		chomp $line;
		if($line ne "" || length($line) != 0){
			if($line =~ /^\s*([#])/){
				next;
			}
			#get the key and the value
			my @params = split("=", $line);
			my $key = $params[0];
			my $val = $params[1];
			if($val eq "" || length($val) == 0){
				croak("Invalid key=value pairs in params.conf file. Key $key has no values.");
			}
			if($key eq 'distribution'){
				if($val eq 'powerlaw'){
					$parameters{"distribution"} = "powerlaw";
				}elsif($val eq 'lognormal'){
					$parameters{"distribution"} = "lognormal";
				}elsif($val eq 'exponential'){
					$parameters{"distribution"} = "exponential";
				}
			}elsif($key eq 'xmin'){
				$parameters{"xmin"} = $val;
			}elsif($key eq 'alpha'){
				$parameters{"alpha"} = $val;
			}elsif($key eq 'mu'){
				$parameters{"mu"} = $val;
			}elsif($key eq 'sigma'){
				$parameters{"sigma"} = $val;
			}elsif($key eq 'lambda'){
				$parameters{"lambda"} = $val;
			}elsif($key eq 'class'){
				$parameters{"class"} = $val;
			}elsif($key eq 'sample'){
				$parameters{"sample"} = $val;
			}elsif($key eq 'read'){
				$parameters{"read"} = $val;
			}elsif($key eq 'evalue'){
				$blast_params{"evalue"} = $val;
			}elsif($key eq 'word'){
				$blast_params{"word_size"} = $val;
			}elsif($key eq 'identity'){
				$blast_params{"perc_identity"} = $val;
			}elsif($key eq 'reward'){
				$blast_params{"reward"} = $val;
			}elsif($key eq 'penalty'){
				$blast_params{"penalty"} = $val;
			}elsif($key eq 'download'){
				$parameters{"download"} = $val;
			}
		}
	}
	$parameters{"blast"} = \%blast_params;
	return %parameters;
}

#This function downloads the viral genome from the NCBI database
sub downloadViralGenomes{
	my($class) = @_;
	`./data.sh`;
}

# Read the genome file names that were downloaded from NCBI
sub readGenomesFileNames{
	my($class) = @_;
	my @genome_file_names = ();
	#read all file names
    unless(open(FileHandle, "data/viruses.txt")){
    	croak("Cannot open file data/viruses.txt. Make sure it is a valid file path. ");
    }	
    my @data = <FileHandle>;
    close(FileHandle);
    #loop through the content and build the content Array
    foreach my $line (@data){
		chomp $line;
		if($line ne "viruses.txt"){
			push(@genome_file_names, $line);	
		}
	}
	return @genome_file_names;		
}
#Reads length of each genome downloaded from NCBI database
sub getAllViralGenomeLength{
	my($class, $file_names_ref) = @_;
	my %all_genome_length;
	my $write_content = "";
	# read all the genome files and get the length of the genome
	for(my $i = 0; $i < @$file_names_ref ; $i++){
		my $genome = "data/".$file_names_ref->[$i];
		# if file exists read the contents of the file
		if(-e $genome){
			unless(open(FileHandle, $genome)){
			  	croak("Cannot open file $genome. \n");
			}	
			my @data = <FileHandle>;
			close(FileHandle);
			my $g_length = 0;
			#loop through the content
			foreach my $line (@data){
				chomp $line;
				if($line =~ /^\s*([>])/){
					next;
				}else{
					$g_length += length($line);
				}
			}
			$all_genome_length{$file_names_ref->[$i]} = $g_length;
			$write_content = $write_content.$file_names_ref->[$i]."\t".$g_length."\n";
		}
	}
	
	# write the genome length into a file
	unless(open(OutFileHandle, ">data/genome_length.txt")){
		croak("Failed to create genome_length.txt file \n");
	}
	unless(print OutFileHandle $write_content){
		croak("Failed to write the content to genome_length.txt file \n");
	}
	
	return %all_genome_length;	
}

#Generates viral community based on the parameter.
sub generateCommunity{
	my($class, $param_ref) = @_;
	my($dist_ref);
	my $distribution = $param_ref->{"distribution"};
	my $community_size = $param_ref->{"class"};
	my $sample_size = $param_ref->{"sample"};
	my $xmin = $param_ref->{"xmin"};
	my $alpha = $param_ref->{"alpha"};
	my $mu = $param_ref->{"mu"};
	my $sigma = $param_ref->{"sigma"};
	my $lambda = $param_ref->{"lambda"};
	
	my (@community);
	if($distribution eq "powerlaw"){
		$dist_ref = Weitz::BioDiversity::DistributionFactory->instantiate("Powerlaw", $community_size);
		@community = $dist_ref->getCommunity(xmin=>$xmin, alpha=>$alpha);
	}elsif($distribution eq "lognormal"){
		$dist_ref = Weitz::BioDiversity::DistributionFactory->instantiate("Lognormal", $community_size);
		@community = $dist_ref->getCommunity(mu=>$mu, sigma=>$sigma);
	}elsif($distribution eq "exponential"){
		$dist_ref = Weitz::BioDiversity::DistributionFactory->instantiate("Exponential", $community_size);
		@community = $dist_ref->getCommunity(lambda=>$lambda);
	}
	return \@community;	
}

# Generates samples from the community
sub generateSample{
	my($class, $sample_size, $comm_aref) = @_;
	my $sample_ref = Weitz::BioDiversity::SampleCommunity->new(community=>$comm_aref, sample_size=>$sample_size);
	my @sample = $sample_ref->getSample();
	my %samples = ("community"=>$comm_aref, "sample"=>\@sample);
	return %samples;	
}

#Simulate reads
sub simulateReads{
	my ($class, $selected_genome_aref, $sample_aref, $read_len) = @_;
	#replacing Metasim with simple read simulator
	#my $read_file = Weitz::BioDiversity::MetaSim->generate454Reads($selected_genome_aref, $sample_aref, $read_len);
	my $read_file = Weitz::BioDiversity::SimpleReadGenerator->generateReads($selected_genome_aref, $sample_aref, $read_len);
	my $count = 0;
	#read the file
    unless(open(FileHandle, $read_file)){
    	croak("Cannot open file $read_file. Make sure it is a valid file path. ");
    }	
    my @data = <FileHandle>;
    close(FileHandle);
    #loop through the content and build the content Array
    foreach my $line (@data){
		chomp $line;
		if($line =~ /^\s*([>])/){
			$count++;	
		}
	}
	my %sim_result = ("read_file"=>$read_file, "read_count"=>$count);
	return %sim_result;	
}

sub makeBlastDB{
	my ($class, $read_file) = @_;
	my $blast_db = Weitz::BioDiversity::BlastPlus->makeBlastDB($read_file);
	return $blast_db;
}

sub runBlast{
	my ($class, $db, $read_file, $sample_size, $blast_param_href) = @_;
	my $blast_out_file = Weitz::BioDiversity::BlastPlus->getBlastMatches($db, $read_file, $sample_size, $blast_param_href);
	return $blast_out_file;
}

sub getBlastResult{
	my ($class, $blast_file) = @_;
	my %results = ();
	my $m_count = 0;
	my $FP = 0;
	unless(open(FileHandle, $blast_file)){
	  	croak("Cannot open file $blast_file. \n");
	}	
	my @data = <FileHandle>;
	close(FileHandle);	
	my $pre_q_read = "";
	my $pre_s_read = "";
	foreach my $line (@data){
		chomp $line;
		#get read names i.e. subject and query
		my $query = (split /\s+/, $line)[0];
		my $subject = (split /\s+/, $line)[1];
		my $q_read = substr($query,0, index($query, '|'));
		my $s_read = substr($subject,0, index($subject, '|'));
		if($q_read eq $s_read || ($pre_q_read eq $q_read && $pre_s_read eq $s_read)){
			$pre_q_read = $q_read;
			$pre_s_read = $s_read;
			next;
		}else{
			$m_count++;
			$pre_q_read = $q_read;
			$pre_s_read = $s_read;			
			my $q_gi = substr($query, rindex($query, '|')+1);
			my $s_gi = substr($subject, rindex($subject, '|')+1);
			if($q_gi ne $s_gi){
				$FP++;
			}
		}
	}
	$results{"TM"} = $m_count;
	$results{"FP"} = $FP;
	return %results;	
}

sub estimateFalseValues{
	my ($class, $blast_href, $sample_aref, $read_count) = @_;
	my $match_count = $blast_href->{"TM"};
	my $false_positive = $blast_href->{"FP"};
	my $poss_match = 0;
	# compute totl possible right matches
	foreach my $item (@$sample_aref){
		if($item != 0){
			$poss_match += $item*($item-1);	
		}
	}
	my $fn_prob = 1 - (($match_count-$false_positive)/$poss_match);
	my $fp_prob = $false_positive/($read_count*($read_count-1));
	$fn_prob = sprintf("%.6f",$fn_prob);
	$fp_prob = sprintf("%.6f",$fp_prob);
	
	my %false_values = ("positive"=>$fp_prob, "negative"=>$fn_prob);
	return %false_values;	
}

# Estimate sample similarity
sub estimateSimilarity{
	my ($class, $community_aref, $match_count, $read_count) = @_;
	my $com_similarity = 0;
	for(my $i = 0; $i < @$community_aref; $i++){
		$com_similarity += $community_aref->[$i]**2;
	}
	# testing with only match_count insteaad of 2*match_count
	my $sample_similarity = ($match_count) /($read_count*($read_count-1));
	$com_similarity = sprintf("%.6f",$com_similarity);
	$sample_similarity = sprintf("%.6f",$sample_similarity);
	my %similarity = ("CS"=>$com_similarity, "SS"=>$sample_similarity);
	
	return %similarity;
}

sub getSortedKeysOfHashByDescValue{
	my ($class, %hash) = @_;
	
	my @skeys = (sort { $hash{$b} <=> $hash{$a} } keys %hash);
	return @skeys;
}

sub spliceArray{
	my (@array, $num) = @_;
	for(my $i = 0; $i < @array;){
		my $s;
		if($i == 0){
			$s = $i;
		}else{
			$s = $i + 1;
		}
		if($i == 0){
			$i = ($i + $num) - 1; 	
		}else{
			$i = $i + $num;
		}
		
		my $e = $i;
		if($e >= @array){$e = @array - 1;}
		if($e >= $s){
			my @fa = @array[$s..$e];
			print "Spliced array @fa \n";
		}
	}	
}

sub cleanRunDir{
	my ($class, $file_path) = @_;
	#get the directory name from the file_path
	my $out = substr($file_path,0, index($file_path, '/'));	
	`rm -r $out`;			
}

1;

=head1 Name
	File: Weitz::BioDiversity::Util.pm
=head1 Description
	This class provides utility methods to run viral diversity estimation
=head1 Synopsis

=head1 Author
	
	Abhiram Das
	
=head1 Copyright

	Copyright (c) 2011 Weitz Lab Georgia Institute of Technology

=cut
