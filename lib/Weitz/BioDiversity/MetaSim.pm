package Weitz::BioDiversity::MetaSim;

use strict;
use warnings;
use Carp;
use Weitz::BioDiversity::Time;

# Generate 454 reads
sub generate454Reads
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
		if($num_reads > 0){
			my $cycles = int($read_length * 0.4);
			#Run metasim with --threads option in multiprocessor machine
			`MetaSim cmd -d $out_dir -4 --454-cycles $cycles -r$num_reads --threads 5 data/$genome_file_name >> $out_dir/simulation.log`;
			
			# rewrite the reads with unique header information
			my $outFile = "$out_dir/".$base_name."-454.fna";
			if(-e $outFile){
				unless(open(FileHandle, $outFile)){
				  	croak("Cannot open file $outFile. \n");
				}	
				my @data = <FileHandle>;
				close(FileHandle);
				#loop through the content and build the content Array
				foreach my $line (@data){
					chomp $line;
					if($line =~ /^\s*([>])/){
						# format the sequence identifier
						# get the gi number from the sequence
						my $index = index($line, 'SOURCES={GI=');
						my $ind_1 = index($line, ',');
						my $sl = length('SOURCES={GI=');
						my $gi = substr($line, $index+$sl, $ind_1-($index+$sl));
						$outString = $outString.">read.".$count."|ref|".$base_name."|GI|".$gi."\n";
						$count++;	
					}else{
						$outString = $outString.$line."\n";
					}
				}
				# Write the reads to a separate file
				`cat $out_dir/$base_name-454.fna >> $out_dir/all-454.fna`;					
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
	File: MetaSim.pm
=head1 Description

=head1 Synopsis

=head1 Author
	
	Ahiram Das
	
=head1 Copyright

	Copyright (c) 2011 Weitz Lab Georgia Institute of Technology

=cut