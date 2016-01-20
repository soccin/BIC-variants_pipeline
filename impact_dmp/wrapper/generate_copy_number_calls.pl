#!/opt/bin/perl
# generate_copy_number_calls.pl ---
# Author: Zehir <zehira@phoenix-h1>
# Created: 04 Feb 2014
# Version: 0.01

use warnings;
use strict;

#my $folder     = $ARGV[0];
#my $sampleList = $ARGV[0];
#my ( %sampleList, %cBioLookup );

#open IN, "<", $sampleList;
#while (<IN>) {
#	chomp;

#	my @line = split("\t");
#	$sampleList{ $line[2] } = 1;
#	$cBioLookup{ $line[2] } = $line[1];

	#print "$line[2]\n";

#}




#my @files = glob("/dmp/hot/zehira/cBioPortalData/copyNumberGeneFiles/*copynumber_segclusp.genes.txt");


die "Error: Please specify input files\n" if(scalar(@ARGV) == 0);
my @files = @ARGV;

my ( %results, @header, %toPrint, %samples );
my $i = 0;


foreach my $file (@files) {

	#print "$file\n";
	#my ($run) = $file =~ /.*\/HiSeqValidationRun(.*)_copynumber.*/;
	#print "$run\n";

	open (IN, $file) or die "Error: Can't open input file $file\n";
	while (<IN>) {
		chomp;
		next if (/sample/);
		my @line   = split("\t");
		my $sample = $line[0];
		#if ( exists( $sampleList{$sample} ) ) {
			#my $sampleID = $cBioLookup{$sample};

			my $gene = $line[3];

			my $significant = $line[9];
			my $FC          = $line[7];
			next if ( $gene =~ /Tiling/ || $gene =~ /FP/ );

			my $call;
			if ( $significant == 1 ) {
				if ( $FC > 0 ) {
					$call = 2;
				}
				elsif ( $FC < 0 ) {
					$call = -2;
				}
			}
			else {
				$call = 0;
			}
			#if ( $i % 342 == 0 ) {
			#	if ( exists( $results{$sample} ) ) {
			#		print STDERR "$sample\n";
			#	}
			#}
			#push @{ $results{$sampleID} }, $gene . ":" . $call;
			push @{ $results{$sample} }, $gene . ":" . $call;
			$i++;
		#}
	}
}

#my (%test);
foreach my $sample ( sort keys %results ) {

	#print "$sample\n";

	push @header, $sample;

	my @results = @{ $results{$sample} };
	my $len = scalar(@results);
	if($len > 341){
		print STDERR "WARNING: This sample is replicated: $sample\t$len !!!!!\n";
	}
	foreach (@results) {
		my ( $gene, $call ) = split(":");
		if ( $call == 2 ) {

		}
		elsif ( $call == -2 ) {

		}
		push @{ $toPrint{$gene} }, $call;
	}

}

print "Hugo_Symbol\t";
my $b = join( "\t", @header );
print "$b\n";

#print scalar(@header), "\t$i\n";

foreach my $gene ( sort keys %toPrint ) {
	print "$gene\t";
	my $a = join( "\t", @{ $toPrint{$gene} } );
	print "$a\n";

}

__END__

=head1 NAME

generate_copy_number_calls.pl - Describe the usage of script briefly

=head1 SYNOPSIS

generate_copy_number_calls.pl [options] args

      -opt --long      Option description

=head1 DESCRIPTION

Stub documentation for generate_copy_number_calls.pl, 

=head1 AUTHOR

Zehir, E<lt>zehira@phoenix-h1E<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014 by Zehir

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut
