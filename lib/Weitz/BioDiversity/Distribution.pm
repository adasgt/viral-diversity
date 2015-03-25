package Weitz::BioDiversity::Distribution;

use strict;
use warnings;
use Carp;

sub new
{
	my ($class, %arg) = @_;
	my $self = bless{
		_num_classes => $arg{num_classes} || croak("Error: Number of classes required."),
	}, $class;
	return $self;
}

# An abstract method that all the inheriting classes will override.
sub getCommunity
{
	my $callers_class = ref($_[0]);
	if($callers_class eq "Weitz::BioDiversity::Distribution")
	{
		croak "Weitz::BioDiversity::Distribution is an abstract base class.";
	}else
	{
		croak "Class $callers_class inherited the abstract base class Weitz::BioDiversity::Distribution 
		but failed to implement the getCommunity method. Attempt to call ${callers_class}::getCommunity";
	}
}

sub get_num_classes {return $_[0]->{_num_classes};}

sub set_num_classes
{
	my ($class, $num_classes) = @_;
	$class->{_num_classes} = $num_classes;
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