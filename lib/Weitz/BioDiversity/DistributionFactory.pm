package Weitz::BioDiversity::DistributionFactory;

use strict;
use warnings;
use Carp;

sub instantiate
{
	my $self = shift;
	my $request_type = shift;
	my $num_classes = shift;
	my $class = "Weitz::BioDiversity::$request_type";
	return $class->new(num_classes=>$num_classes);
}

1;

=head1 Name

	Weitz::BioDiversity::DistributionFactory
	
=head1 Description


=head1 Synopsis


=head1 Author

	Abhiram Das

=head1 Copyright

	Copyright (c) 2011 Weitz Lab Georgia Institute of Technology
	
=cut