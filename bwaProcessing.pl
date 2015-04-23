#!/usr/bin/perl

use strict;

### input is samse/sampe file from bwa

### NOTE: WILL ONLY KEEP READS THAT ARE CONSIDERED PROPER PAIR BASED ON BWA BIT FLAG
### keeps only unique single end reads

### this will not keep singletons that have mapQ == 0; 

my $file = shift;
open(IN, "$file");
open(OUT, ">$file\_filtered.sam");
my %pairs = ();
while(my $line = <IN>){
    chomp $line;


    if($line =~ /^\@/){
	if($line !~ /\@HD/){
	    print OUT "$line\n";
	}
	next;
    }

    my @data = split(/\s+/, $line);

    my $readPairedInSequencing = $data[1] & 1;
    my $temp2 = $data[1]>>1;
    my $readMappedInProperPair = $temp2 & 1;
    my $temp4 = $data[1]>>2;
    my $queryUnmapped = $temp4 & 1;
    my $temp8 = $data[1]>>3;
    my $mateUnmapped = $temp8 & 1;
    my $temp16 = $data[1]>>4;
    my $strand = $temp16 & 1;
    #my $temp64 = $data[1]>>6;
    #my $firstRead = $temp64 & 1;

    if($queryUnmapped){
	next;
    }

    my $flag = 0;
    if($readPairedInSequencing){
	if($mateUnmapped){
	    ### if mate unmapped, keep read if unique alignment
	    ### for some reason, bwa doesn't always give correct pair info so must change it here
	    ### for some reason, XT:A:U aren't always unique alignments so using mapping qual also
	    if($line =~ /XT:A:U/ && $data[4] > 0){
		my $singlified = singled(@data);		
		print OUT "$singlified\n";
		next;
	    }
	}
	else{
	    ### if unique, keep end
	    if($line =~ /XT:A:U/ && $data[4] > 0){
		$flag = 1;
	    }
	    else{
		if($readMappedInProperPair){
		    ### if repetitive, keep only if proper pair
		    $flag = 1;
		}
	    }
	}
    }
    else{
	### if single end, keep if unique alignment
	if($line =~ /XT:A:U/ && $data[4] > 0){
	    print OUT "$line\n";
	    next;
	}
    }

    if($flag == 1){
	if(!$pairs{$data[0]}){
	    my @temp = ();
	    push @temp, $line;
	    $pairs{$data[0]} = [@temp];
	}
	else{
	    my @temp2 = @{$pairs{$data[0]}};
	    push @temp2, $line;
	    $pairs{$data[0]} = [@temp2];
	}
    }
}
close IN;

foreach my $key (keys %pairs){
    if(scalar(@{$pairs{$key}}) == 1){
	### removes repetitive single ends
	my @data = split(/\s+/, @{$pairs{$key}}[0]);
	if(@{$pairs{$key}}[0] =~ /XT:A:U/ && $data[4] > 0){
	    my @data = split(/\s+/, @{$pairs{$key}}[0]);
	    my $readPairedInSequencing = $data[1] & 1;
	    if($readPairedInSequencing){
		### paired reads whose mates were removed by markDups
		### or uniques with the other end not placing in proper pair
		### change read alignment information to that of a single end
				
		my $singlified = singled(@data);		
		print OUT "$singlified\n";
	    }
	}
    }
    elsif(scalar(@{$pairs{$key}}) == 2){
	my @first = split(/\s+/, @{$pairs{$key}}[0]);
	my @second = split(/\s+/, @{$pairs{$key}}[1]);

	my $temp21 = $first[1]>>1;
	my $readMappedInProperPair1 = $temp21 & 1;

	my $temp22 = $second[1]>>1;
	my $readMappedInProperPair2 = $temp22 & 1;

	### if both ends are unique, then only keep if on same chr, correct insert length, and in proper orientation;
	if((@{$pairs{$key}}[0] =~ /XT:A:U/ && $first[4] > 0) && (@{$pairs{$key}}[1] =~ /XT:A:U/ && $second[4] > 0)){
	    if($readMappedInProperPair1 && $readMappedInProperPair2){
		print OUT "@{$pairs{$key}}[0]\n";
		print OUT "@{$pairs{$key}}[1]\n";
	    }
	}
	elsif((@{$pairs{$key}}[0] =~ /XT:A:U/ && $first[4] > 0) && (@{$pairs{$key}}[1] !~ /XT:A:U/ || $second[4] != 0)){
	    ### if one end is unique and other isn't
	    ### keep repetitive end if proper pair
	    if($readMappedInProperPair1 && $readMappedInProperPair2){
		print OUT "@{$pairs{$key}}[0]\n";
		print OUT "@{$pairs{$key}}[1]\n";
	    }
	    else{
		my $singlified = singled(@first);		
		print OUT "$singlified\n";
	    }
	}
	elsif((@{$pairs{$key}}[0] !~ /XT:A:U/ || $first[4] != 0) && (@{$pairs{$key}}[1] =~ /XT:A:U/ && $second[4] > 0)){
	    if($readMappedInProperPair1 && $readMappedInProperPair2){
		print OUT "@{$pairs{$key}}[0]\n";
		print OUT "@{$pairs{$key}}[1]\n";
	    }
	    else{
		my $singlified = singled(@second);		
		print OUT "$singlified\n";
	    }
	}
    }
}
close OUT;  


sub singled {
    my @alignment = @_;

    my $temps = $alignment[1]>>4;
    my $strand = $temps & 1;

    if($strand){
	$alignment[1] = 16;
    }
    else{
	$alignment[1] = 0;
    }
    
    $alignment[6] = "*";
    $alignment[7] = 0;
    $alignment[8] = 0;

    my $newLine = join("\t", @alignment);
    return $newLine;
}
