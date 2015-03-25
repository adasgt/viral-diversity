package Weitz::BioDiversity::BlastPlus;

use strict;
use warnings;
use Carp;
#use Bio::Seq;
#use Bio::SeqIO;
#use Bio::Tools::Run::StandAloneBlastPlus;
#use Bio::Search::Result::BlastResult;

#
# Uses BioPerl BlastPlus module to blast the output
# 
sub getPairWiseMatches
{
	my($class, $seq_file, $param_ref) = @_;
	
	my $match_count = 0;
	
	my $fac = Bio::Tools::Run::StandAloneBlastPlus->new(-create => 1); # -db_name=> 'virusdb', -db_data=>$seq_file,
	my $seq_obj = Bio::SeqIO->new(-file => $seq_file, -format => "fasta", -alphabet => "dna");
	
	#loop through each seq in the seq-obj and blast it with the next sequence 
	my @seq_list;
	while(my $seq = $seq_obj->next_seq()){
		push(@seq_list, $seq);
	}
	# read blast paramas
	my $eval = $param_ref->{"evalue"} + 0;
	my $word = $param_ref->{"word_size"} + 0;
	my $identity = $param_ref->{"perc_identity"} + 0;
	my $reward = $param_ref->{"reward"} + 0;
	my $penalty = $param_ref->{"penalty"} + 0;
       
	for(my $i = 0; $i < @seq_list; $i++){
	    for(my $j = $i+1; $j < @seq_list; $j++){
			my $result = $fac->bl2seq(-method=>"blastn", -query => $seq_list[$j], -subject => $seq_list[$i], -method_args => ['-evalue' => $eval, '-num_threads' =>2, '-word_size'=>$word, '-perc_identity'=>$identity, '-penalty'=>$penalty, '-reward'=> $reward]);
			$fac->rewind_results();
			if($result = $fac->next_result){
				if(my $hit = $result->next_hit){
					if (my $hsp = $hit->next_hsp){
						$match_count++;
					}
				}
			}
		}
	}

	$fac->cleanup();
	return $match_count;
}

sub getBlastMatches{
	my($class, $blast_db, $seq_file, $num_seq, $param_ref) = @_;
	
	# read blast paramas
	my $eval = $param_ref->{"evalue"} + 0;
	my $word = $param_ref->{"word_size"} + 0;
	my $identity = $param_ref->{"perc_identity"} + 0;
	my $reward = $param_ref->{"reward"} + 0;
	my $penalty = $param_ref->{"penalty"} + 0;	
	#get the directory name from the seq_file
	my $out = substr($blast_db,0, index($blast_db, '/'));	
	
	my $blastout = "$out/blast_out.out";
	
	#run blast
	`blastn -db $blast_db -query $seq_file -evalue $eval -word_size $word -perc_identity $identity -penalty $penalty -reward $reward -out $blastout -outfmt 6 -num_threads 5`;
	return $blastout;
}

sub makeBlastDB{
	my($class, $seq_file) = @_;
	#get the directory name from the seq_file
	my $out = substr($seq_file,0, index($seq_file, '/'));
	my $blast_db = "$out/all_blast.db";
	
	# Build the blast database
	`makeblastdb -in $seq_file -dbtype nucl -out $blast_db -logfile $out/blast_makedb.log`;
	
	return $blast_db;
}

1;

=head1 Name
	File: Weitz::BioDiversity::BlastPlus.pm
=head1 Description

=head1 Synopsis

=head1 Author
	
	Abhiram Das
	
=head1 Copyright

	Copyright (c) 2011 Weitz Lab Georgia Institute of Technology

=cut
