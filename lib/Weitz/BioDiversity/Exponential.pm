package Weitz::BioDiversity::Exponential;

use base (Weitz::BioDiversity::Distribution);
use strict;
use warnings;
use Carp;

# Get community for a Exponential distribution
# Formula for Exponential distribution is as follows::
# x = -(ln x)/lambda
sub getCommunity
{
	my ($class, %arg) = @_;
	if((keys %arg) < 1) {croak("Invalid number of arguments passed to Exponential::getCommunity() method.");}
	my $lambda = $arg{$_[1]};
	
	my @community = ();
	#populate the community
	for(my $i = 0; $i < $class->get_num_classes(); $i++){
		$community[$i] = sprintf("%.10f",rand());
	}
	for(my $i = 0; $i < @community; $i++){
		$community[$i] = -log($community[$i]) / $lambda;
	}	
	# Take the cumulative sum of the community
	my $cum_sum = 0;
	for(@community){$cum_sum += $_}
	# Normalize the community by cumulatve sum
	for(my $i = 0; $i < @community; $i++){
		$community[$i] = $community[$i]/$cum_sum;
	} 
	return @community;
}

1;

=head1 Name

	Weitz::BioDiversity::Exponential
	
=head1 Description

	It's a Exponential distribution class implementation. It allows to get the community for a Exponential distribution.

=head1 Synopsis


=head1 Author

	Abhiram Das

=head1 Copyright

	Copyright (c) 2011 Weitz Lab Georgia Institute of Technology
	
=cut