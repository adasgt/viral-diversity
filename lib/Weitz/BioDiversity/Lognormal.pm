package Weitz::BioDiversity::Lognormal;

use base (Weitz::BioDiversity::Distribution);
use strict;
use warnings;
use Carp;
use Math::Random;

# Get community for a Lognormal distribution
# Formula for lognormal distribution is as follows::
# f(x:mu,sigma) = (1/(x*sigma*squareroot(2pi))) exp (-((ln x -mu)^2)/2sigma^2)
sub getCommunity
{
	my ($class, %arg) = @_;
	if((keys %arg) < 2) {croak("Invalid number of arguments passed to Lognormal::getCommunity() method.");}
	my $mu = $arg{$_[1]};
	my $sigma = $arg{$_[3]};
	my @community;
	# replace the following code. Now uses a normal random generator
	for(my $i = 0; $i < $class->get_num_classes(); $i++){
		$community[$i] = sprintf("%.10f",random_normal());
	}
	for(my $i = 0; $i < @community; $i++){
		$community[$i] = exp($mu + ($sigma * $community[$i]));
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

	Weitz::BioDiversity::Lognormal
	
=head1 Description

	It's a Lognormal distribution class implementation. It allows to get the community for a Lognormal distribution.

=head1 Synopsis


=head1 Author

	Abhiram Das

=head1 Copyright

	Copyright (c) 2011 Weitz Lab Georgia Institute of Technology
	
=cut