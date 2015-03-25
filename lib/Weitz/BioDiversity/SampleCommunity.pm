package Weitz::BioDiversity::SampleCommunity;

use strict;
use warnings;
use Carp;
use List::MoreUtils qw 'lastidx';

sub new
{
	my ($class, %arg) = @_;
	my $self = bless{
		_community => $arg{community} || croak("Error: Community is Required."),
		_sample_size => $arg{sample_size} || croak("Error: Sample Size is Required."),
	}, $class;
	return $self;
}

# Get the sample from the community
sub getSample
{
	my($class) = @_;
	# generate cumulative probability distribution of the communiuty
	my @com_dist = ();
	my $com_ref = $class->{_community};
	for(my $i = 0; $i < @$com_ref ; $i++){
		if($i==0){
			$com_dist[$i] = $class->{_community}[$i];
		}else{
			$com_dist[$i] = $class->{_community}[$i] + $com_dist[$i-1];
		}		
	}
	my @temp_sample = ();
	# Create random sample of size _sample_size
	for(my $i = 0; $i < $class->{_sample_size}; $i++){
		$temp_sample[$i] = sprintf("%.10f",rand());
	}
	@temp_sample = sort(@temp_sample);
	my @sample = ();
	# loop through the com_dist array and find out how many elements in the temp_sample fall below each element 
	# of com_dis
	my $pre_index = 0; 	
	for(my $i = 0; $i < @com_dist; $i++){
		my $count = 0;
		my $last_index = 0;
		if($i == 0){
			$last_index = lastidx{$_ <= $com_dist[$i]} @temp_sample;
			$count = ($last_index - $pre_index) +1;
			$pre_index = $last_index;
		}else{
			$last_index = lastidx{$_ <= $com_dist[$i]} @temp_sample;
			$count = ($last_index - $pre_index);
			$pre_index = $last_index;			
		}
		$sample[$i] = $count;
	}
	return @sample;
}


1;

=head1 Name

=head1 Description

=head1 Synopsis

=head1 Author
	
	Ahiram Das
	
=head1 Copyright

	Copyright (c) 2011 Weitz Lab Georgia Institute of Technology

=cut