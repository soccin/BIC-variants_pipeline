#!/opt/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin); 

### takes in picard metrics (tested on hs/is/as/md
### through -metrics
### and/or can give it a file of files (-files)

my $pre = 'TEMP';
my @metrics = ();
my $files = '';
GetOptions ('pre=s' => \$pre,
	    'files=s' => \$files,
	    'metrics|m=s' => \@metrics) or exit(1);

if($files){
    open(FI, "$files") or die "Can't open file $files $!";
    while(<FI>){
	chomp;

	push @metrics, $_;
    }
    close FI;
}
    
my $hflag = 0;

foreach my $met (@metrics){
    if(!-s $met){
	die "METRICS FILE $met DOES NOT EXIST OR IS EMPTY";
    }
}

foreach my $met (@metrics){
    open(MET, $met);
    my $mflag = 0;
    while(my $line=<MET>){
	chomp $line;

	if($line =~ /^\#\# METRICS CLASS/){
	    if($mflag == 0){
		$mflag = 1;
	    }

	    my $header = <MET>;
	    if($hflag == 0){
		print "$header";
		$hflag = 1;
	    }
	}
	else{
	    if($mflag == 0){
		next;
	    }
	    else{
		if($line){
		    print "$line\n";
		}
		else{
		    last;
		}
	    }
	}
    }
    close MET;
}
