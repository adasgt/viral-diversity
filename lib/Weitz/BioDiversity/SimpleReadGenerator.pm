package Weitz::BioDiversity::SimpleReadGenerator;

use strict;
use warnings;
use Carp;
use Weitz::BioDiversity::Time;

# Generate reads
# Arguments:
# genome_name_ref - List of genome files
# sample_ref - List containing number of reads to be generated for each corresponding genome in the genome list
# read_length - Average read length

sub generateReads
{
	my($class, $genome_name_ref, $sample_ref, $read_length) = @_;
	
	my $outString = "";
	my $count = 1;
	my $out_dir = "out_".Weitz::BioDiversity::Time->new()->toShortString();
	#make Directory
	`mkdir $out_dir`;
	# delete files from out directory
	#`rm out/*.*`;
	# get the genome sequence file and generate required number of reads
	for(my $i = 0; $i < @$genome_name_ref; $i++){
		my $genome_file_name = $genome_name_ref->[$i];
		my $base_name = substr($genome_file_name, 0, index($genome_file_name, '.'));
		my $num_reads = $sample_ref->[$i];
		my $seq_head = "";
		my $sequence = "";
		my $seq_length = 0;
		# Read the genome file. Get the header and the genome string. Get the genome length
		# if file exists read the contents of the file
		if(-e "data/".$genome_file_name){
			unless(open(FileHandle, "data/$genome_file_name")){
			  	croak("Cannot open file data/$genome_file_name. \n");
			}	
			my @data = <FileHandle>;
			close(FileHandle);
			#loop through the content
			foreach my $line (@data){
				chomp $line;
				if($line =~ /^\s*([>])/){
					$seq_head = $line;
					next;
				}else{
					$sequence = $sequence.$line;
					$seq_length += length($line);
				}
			}
		}
		
		if($num_reads > 0){
			# Generate num_read random numbers (index in the seq) between 0 and seq_length
			my @seq_index;
			for(my $i = 0; $i < $num_reads; $i++){
				$seq_index[$i] = int rand($seq_length);
			}
			# create suffix for the read header i.e. |GI|#|ref|ref_name
			my $se_gi = substr($seq_head, index($seq_head, "gi")+3, index($seq_head, "ref")-1 - (index($seq_head, "gi")+3));
			#generate reads
			for(my $i=0; $i < $num_reads; $i++){
				# check if seq_index + read_lenghth < sequence length i.e. its not end of the sequence
				# if end of sequence read in reverse direction
				my $read = ">read.".$count."|gi|".$se_gi."|ref|".$base_name."\n";
				if($seq_index[$i] + $read_length < $seq_length){
					$read = $read.substr($sequence, $seq_index[$i], $read_length)."\n";
				}else{
					$read = $read.scalar reverse(substr($sequence, $seq_index[$i] - $read_length, $read_length))."\n";
				}
				#build the output string
				$outString = $outString.$read;
				$count++;
			}
		}
	}
	
	# Write the output string to the file
	unless( open( FileHandle, ">$out_dir/all_reads.fasta") ) {
		croak("Cannot write to file $out_dir/all_reads.fasta ");
	}
	unless(print FileHandle $outString) {
		croak("Cannot write to file $out_dir/all_reads.fasta");
	}
	close(FileHandle);		 
	return "$out_dir/all_reads.fasta";
}


1;

=head1 Name
	File: SimpleReadGenerator.pm
=head1 Description
	Generate simple reads for a given genome
=head1 Synopsis

=head1 Author
	
	Abhiram Das
	
=head1 Copyright

	Copyright (c) 2011 Weitz Lab Georgia Institute of Technology

=cut
