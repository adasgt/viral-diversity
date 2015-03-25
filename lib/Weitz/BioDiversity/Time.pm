package Weitz::BioDiversity::Time;

use strict;
use warnings;
use Carp;

my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);

sub new{
	my ($class) = @_;
	my $self = bless{}, $class;
	($self -> {second}, $self -> {minute}, $self -> {hour}, $self -> {dayOfMonth}, 
	$self -> {month}, $self -> {yearOffset}, $self -> {dayOfWeek}, $self -> {dayOfYear},$self -> {dayLightSavings}) = localtime(time);
	
	$self -> {year} = 1900 + $self -> {yearOffset};
	
	return $self;
}

sub getLocalTime{
	my ($self) = @_;
	
	my $currentTime = $self -> {hour}.":".$self -> {minute}.":".$self -> {second}.", ". $weekDays[$self -> {dayOfWeek}]." ".
						$months[$self -> {month}]." ".$self -> {dayOfMonth}.", ".$self -> {year};
	return $currentTime;
}

sub getGMTime{
	my ($self) = @_;
	($self -> {second}, $self -> {minute}, $self -> {hour}, $self -> {dayOfMonth}, 
	$self -> {month}, $self -> {yearOffset}, $self -> {dayOfWeek}, $self -> {dayOfYear},$self -> {dayLightSavings}) = gmtime(time);
	
	$self -> {year} = 1900 + $self -> {yearOffset};
	my $currentTime = $self -> {hour}.":".$self -> {minute}.":".$self -> {second}.", ". $weekDays[$self -> {dayOfWeek}]." ".
				$months[$self -> {month}]." ".$self -> {dayOfMonth}.", ".$self -> {year};
	return $currentTime;
}

sub getTime{
	my ($self) = @_;
	
	my $currentTime = ($self->{hour}*60*60 + $self->{minute}*60 + $self->{second}) * 1000;
	return $currentTime;	
}

sub toString{
	my ($self) = @_;
	my $currentTime = $months[$self -> {month}]."_".$self -> {dayOfMonth}."_".$self -> {year}."_h".$self -> {hour}."_m".$self -> {minute}."_s".$self -> {second};
	return $currentTime;
}

sub toShortString{
	my ($self) = @_;
	my $currentTime = $self -> {month}.$self -> {dayOfMonth}.$self -> {year}.$self -> {hour}.$self -> {minute}.$self -> {second};
	return $currentTime;
}


1;


=head1 Name

=head1 Description

=head1 Synopsis

=head1 Author
	
	Abhiram Das
	
=head1 Copyright

	Copyright (c) 2011 Weitz Lab Georgia Institute of Technology

=cut