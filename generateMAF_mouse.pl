#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin); 

### INPUT: vcf file and list of normal/tumor sample pairing information
### OUTPUT: 2 maf files; 

### NOTE: CAN'T RUN ON NODE BECAUSE NO NETWORK ACCESS FOR ONCOTATOR
### NOTE2: DESIGNED TO WORK FOR ONCOTATOR ANNOTATED MAFS

my ($vcf, $pairing, $species, $config, $caller, $normal_sample, $tumor_sample, $delete_temp);
GetOptions ('vcf=s' => \$vcf,
	    'species=s' => \$species,
	    'config=s' => \$config,
	    'caller=s' => \$caller,
	    'normal_sample=s' => \$normal_sample,
	    'tumor_sample=s' => \$tumor_sample,
	    'delete_temp' => \$delete_temp,
	    'pairing=s' => \$pairing);


my $somatic = 'UNPAIRED';
if(!$vcf){
    die "YOU MUST PROVIDE A VCF\n";
}

if(!-e $vcf){
    die "$vcf DOES NOT EXIST";
}

if($pairing){
    $somatic = 'PAIRED';
    if(!-e $pairing){
	die "$pairing DOES NOT EXIST";
    }
}

if($species !~ /mm9|mm10/i){
    die "Only support for mm9 or mm10";
}

if($caller !~ /unifiedgenotyper|ug|haplotypecaller|hc|mutect|varscan|somaticsniper/i){
    die "Only support for unifiedgenotyper(ug), haplotypecaller(hc), mutect, varscan, and somaticsniper";
}

my $PYTHON = '';
open(CONFIG, "$config") or warn "CAN'T OPEN CONFIG FILE $config SO USING DEFAULT SETTINGS";
while(<CONFIG>){
    chomp;

    my @conf = split(/\s+/, $_);
    if($conf[0] =~ /python/i){
	if(!-e "$conf[1]/python"){
	    die "CAN'T FIND python IN $conf[1] $!";
	}
	$PYTHON = $conf[1];
    }
}
close CONFIG;


print "converting to MAF and filtering snp calls for coverage/qual\n";
if($caller =~ /unifiedgenotyper|ug|haplotypecaller|hc/i){
    if($pairing){
	print "$PYTHON/python $Bin/maf/vcf2maf0.py -i $vcf -c $caller -o $vcf\_PAIRED.maf -p $pairing\n";
	print "cat $vcf\_PAIRED.maf | $PYTHON/python $Bin/maf/pA_qSomHC.py | /usr/bin/tee $vcf\_$somatic\_maf0.txt\n\n";

	`$PYTHON/python $Bin/maf/vcf2maf0.py -i $vcf -c $caller -o $vcf\_PAIRED.maf -p $pairing`;
	`cat $vcf\_PAIRED.maf | $PYTHON/python $Bin/maf/pA_qSomHC.py | /usr/bin/tee $vcf\_$somatic\_maf0.txt > $vcf\_$somatic\_maf0.log 2>&1`;
    }
    else{
 	print "$PYTHON/python $Bin/maf/vcf2maf0.py -i $vcf -c $caller -o $vcf\_UNPAIRED.maf\n";
	print "cat $vcf\_UNPAIRED.maf | $PYTHON/python $Bin/maf/pA_Cov+noLow.py | /usr/bin/tee $vcf\_$somatic\_maf0.txt\n\n";

	`$PYTHON/python $Bin/maf/vcf2maf0.py -i $vcf -c $caller -o $vcf\_UNPAIRED.maf`;
	`cat $vcf\_UNPAIRED.maf | $PYTHON/python $Bin/maf/pA_Cov+noLow.py | /usr/bin/tee $vcf\_$somatic\_maf0.txt > $vcf\_$somatic\_maf0.log 2>&1`;
   }
}
else{
    my $af = '';
    if($caller =~ /mutect|mu/i){
	my $mutext = $vcf;
	$mutext =~ s/\.vcf$/\.txt/;
	$af = "-aF $mutext";
    }

    my $n_sample = '';
    my $t_sample = '';
    if($normal_sample && $tumor_sample){
	$n_sample = "-n $normal_sample";
	$t_sample = "-t $tumor_sample";
    }

    print "$PYTHON/python $Bin/maf/vcf2maf0.py -i $vcf $af -c $caller -o $vcf\_$somatic\_maf0.txt -p $pairing $n_sample $t_sample\n\n";
    `$PYTHON/python $Bin/maf/vcf2maf0.py -i $vcf $af -c $caller -o $vcf\_$somatic\_maf0.txt -p $pairing $n_sample $t_sample`;

    if($caller =~ /mutect|mu/i){
	`/bin/mv $vcf\_$somatic\_maf0.txt $vcf\_$somatic\_UNFILTERED.txt`;
	print "$PYTHON/python $Bin/rescue/DMP_rescue.py <$vcf\_$somatic\_UNFILTERED.txt> $vcf\_$somatic\_rescued.txt 2> $vcf\_$somatic\_maf0_rescue.log";
	`$PYTHON/python $Bin/rescue/DMP_rescue.py <$vcf\_$somatic\_UNFILTERED.txt> $vcf\_$somatic\_rescued.txt 2> $vcf\_$somatic\_maf0_rescue.log`;
	
	open(MUT, "$vcf\_$somatic\_rescued.txt");
	open(OUT, ">$vcf\_$somatic\_maf0.txt");
	my $header = <MUT>;
	print OUT "$header";
	my @columns = split(/\t/, $header);
	my $index_keep = -1;
	
	for(my $i=0; $i<scalar(@columns); $i++){
	    if($columns[$i] eq 'MUT_KEEP'){
		$index_keep = $i;
	    }
	}
	
	if($index_keep == -1){
	    die "Can find MUT_KEEP column in $vcf\_$somatic\_UNFILTERED.txt";
	}
	
	while(my $line=<MUT>){
	    chomp $line;
	    
	    my @data = split(/\t/, $line);
	    if($data[$index_keep] eq 'KEEP'){
		print OUT "$line\n";
	    }
	}
	close MUT;
	close OUT;
    }
}

print "filtering for snpeff functioal annotations\n";
print "cat $vcf\_$somatic\_maf0.txt | $PYTHON/python $Bin/maf/pA_Functional_Snpeff.py > $vcf\_$somatic\_maf1.txt\n\n";
`cat $vcf\_$somatic\_maf0.txt | $PYTHON/python $Bin/maf/pA_Functional_Snpeff.py > $vcf\_$somatic\_maf1.txt`;

print "converting to new MAF format using species $species ... \n";
print "$PYTHON/python $Bin/maf/oldMAF2tcgaMAF.py $species $vcf\_$somatic\_maf1.txt $vcf\_$somatic\_TCGA_MAF.txt\n\n";
### Convert old (NDS) MAF to official TCGA MAF
### NOTE; DON'T FORGET <> AROUND INPUT FILE (maf0.txt)
`$PYTHON/python $Bin/maf/oldMAF2tcgaMAF.py $species $vcf\_$somatic\_maf1.txt $vcf\_$somatic\_TCGA_MAF.txt`;

if($delete_temp){
    `/bin/rm $vcf\_PAIRED.maf $vcf\_$somatic\_maf0.txt $vcf\_$somatic\_maf0.log $vcf\_UNPAIRED.maf $vcf\_$somatic\_UNFILTERED.txt $vcf\_$somatic\_rescued.txt $vcf\_$somatic\_maf1.txt`;
}
