#!/opt/bin/perl

package Cluster;
use strict;
use Data::Dumper;

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self  = {};

    $self->{submit} = undef;
    $self->{job_name} = undef;
    $self->{job_hold} = undef;
    $self->{cpu} = undef;
    $self->{mem} = undef;
    $self->{internet} = undef;
    $self->{cluster_out} = undef;
    $self->{cluster_error} = undef;

    bless $self, $class;
    return $self;
}


sub submit {
    my $self = shift;
    if (@_){
	$self->{submit} = shift;
    }
    return $self->{submit};
}

sub job_name {
    my $self = shift;
    if (@_){
	$self->{job_name} = shift;
    }
    return $self->{job_name};
}

sub job_hold {
    my $self = shift;
    if (@_){
	$self->{job_hold} = shift;
    }
    return $self->{job_hold};
}

sub cpu {
    my $self = shift;
    if (@_){
	$self->{cpu} = shift;
    }
    return $self->{cpu};
}

sub mem {
    my $self = shift;
    if (@_){
	$self->{mem} = shift;
    }
    return $self->{mem};
}

sub internet {
    my $self = shift;
    if (@_){
	$self->{internet} = shift;
    }
    return $self->{internet};
}

sub cluster_out {
    my $self = shift;
    if (@_){
	$self->{cluster_out} = shift;
    }
    return $self->{cluster_out};
}

sub cluster_error {
    my $self = shift;
    if (@_){
	$self->{cluster_error} = shift;
    }
    return $self->{cluster_error};
}

1;
