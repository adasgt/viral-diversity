package Weitz::BioDiversity::Powerlaw;

use base (Weitz::BioDiversity::Distribution);
use strict;
use warnings;
use Carp;

# Get community for a powerlaw distribution
sub getCommunity
{
	my ($class, %arg) = @_;
	if((keys %arg) < 2) {croak("Invalid number of arguments passed to Powerlaw::getCommunity() method.");}
	my $xmin = $arg{$_[1]};
	my $alpha = $arg{$_[3]};
	
	my @community = ();
	#populate the community
	for(my $i = 0; $i < $class->get_num_classes(); $i++){
		$community[$i] = sprintf("%.10f",rand());
	}
	for(my $i = 0; $i < @community; $i++){
		$community[$i] = ((1-$community[$i]) ** (-1 / ($alpha-1))) * $xmin;
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

	Weitz::BioDiversity::Powerlaw
	
=head1 Description

	It's a Powerlaw distribution class implementation. It allows to get the community for a power law distribution.

=head1 Synopsis


=head1 Author

	Abhiram Das

=head1 Copyright

	Copyright (c) 2011 Weitz Lab Georgia Institute of Technology
	
=cut