#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin); 
use File::Basename;

### INPUT: vcf file and list of normal/tumor sample pairing information
### OUTPUT: 2 maf files; 

## CONSTANT FOR VEP
my $VEP_COLUMN_NAMES = "Center,Verification_Status,Validation_Status,Mutation_Status,Sequencing_Phase,Sequence_Source,Validation_Method,Score,BAM_file,Sequencer,Tumor_Sample_UUID,Match_Norm_Sample_UUID,Caller";

my ($vcf, $pairing, $patient, $bam_dir, $species, $config, $caller, $normal_sample, $tumor_sample, $delete_temp);
GetOptions ('vcf=s' => \$vcf,
	    'species=s' => \$species,
	    'config=s' => \$config,
	    'caller=s' => \$caller,
            'align_dir=s' => \$bam_dir,
            'patient=s' => \$patient,
	    'normal_sample=s' => \$normal_sample,
	    'tumor_sample=s' => \$tumor_sample,
	    'delete_temp' => \$delete_temp,
	    'pairing=s' => \$pairing) or die;


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

if($patient) {
    if(!-e $patient){
        die "$patient DOES NOT EXIST";
    }
    if(!$bam_dir){
        die "If patient file is given, you must supply alignment directory for fillout.";
    }
}

if($bam_dir){
    if(!-e $bam_dir){
        die "$bam_dir DOES NOT EXIST";
    }
}

if($species !~ /hg19/){ ##    UNCOMMENT LATER   |mm10|mouse/i){
    print "THIS WILL ONLY PRINT OUT A TCGA MAF, NO ANNOATION\n";
}

if($caller !~ /unifiedgenotyper|ug|haplotypecaller|hc|mutect|varscan|somaticsniper/i){
    die "Only support for unifiedgenotyper(ug), haplotypecaller(hc), mutect, varscan, and somaticsniper";
}

my $REF_FASTA = '';
my $HG19_FASTA = '';
my $PYTHON = '';
my $PERL = '';
my $VEP = '';
my $MM10_FASTA = '';
open(CONFIG, "$config") or warn "CAN'T OPEN CONFIG FILE $config SO USING DEFAULT SETTINGS";
while(<CONFIG>){
    chomp;

    my @conf = split(/\s+/, $_);
    if($conf[0] =~ /hg19_fasta/i){
        if(!-e "$conf[1]"){
          die "CAN'T FIND $conf[1] $!";
        }
        $HG19_FASTA = $conf[1];
    }
    elsif($conf[0] =~ /mm10_fasta/i){
        if(!-e "$conf[1]"){
            if($species =~ /^mm10$/i){
                die "CAN'T FIND $conf[1] $!";
            }
        }
        $MM10_FASTA = $conf[1];
    }
    elsif($conf[0] =~ /python/i){
	if(!-e "$conf[1]/python"){
	    die "CAN'T FIND python IN $conf[1] $!";
	}
	$PYTHON = $conf[1];
    }
    elsif($conf[0] =~ /perl/i){
	if(!-e "$conf[1]/perl"){
	    die "CAN'T FIND perl IN $conf[1] $!";
	}
	$PERL = $conf[1];
    }
    elsif($conf[0] =~ /vep/i){
        if(!-e "$conf[1]/variant_effect_predictor.pl"){
            die "CAN'T FIND VEP IN $conf[1] $!";
        }
        $VEP = $conf[1];
    }
    elsif($conf[0] =~ /samtools/i){
        if(!-e "$conf[1]/samtools"){
            die "CAN'T FIND samtools IN $conf[1] $!";
        }
        my $path_tmp = $ENV{'PATH'};
        $ENV{'PATH'} = "$conf[1]:$path_tmp";
    }
    elsif($conf[0] =~ /tabix/i){
        if(!-e "$conf[1]/tabix"){
            die "CAN'T FIND tabix IN $conf[1] $!";
        }
        my $path_tmp = $ENV{'PATH'};
        $ENV{'PATH'} = "$conf[1]:$path_tmp";
    }

}
close CONFIG;

my $NCBI_BUILD = '';
my $VEP_SPECIES = '';
if($species =~ /hg19/i){
    $species = 'hg19';
    $REF_FASTA = "$HG19_FASTA";
    $NCBI_BUILD = "GRCh37";
    $VEP_SPECIES = "homo_sapiens";
}
elsif($species =~ /mouse|^mm10$/i){
    $species = 'mm10';
    $REF_FASTA = "$MM10_FASTA";
    $NCBI_BUILD = "GRCm38";
    $VEP_SPECIES = "mus_musculus";
}

## RIGHT NOW ONLY HG19 IS BEING USED IN THIS SCRIPT

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

print "converting to new MAF format using species $species ... \n";
print "$PYTHON/python $Bin/maf/oldMAF2tcgaMAF.py $species \<$vcf\_$somatic\_maf0.txt\> $vcf\_$somatic\_maf1.txt\n\n";
### Convert old (NDS) MAF to official TCGA MAF
### NOTE; DON'T FORGET <> AROUND INPUT FILE (maf0.txt)
`$PYTHON/python $Bin/maf/oldMAF2tcgaMAF.py $species $vcf\_$somatic\_maf0.txt $vcf\_$somatic\_maf1.txt`;


if($species !~ /hg19/i) { ###   |mm10|mouse/i) { uncomment later!
    print "End of species ambiguous portion of the script.\n";
    exit 0;
}

print "\n#######\n#######\nStarting VEP. \n";
# these are names needed for the "retain-cols" option in VEP
# my $VEP_COLUMN_NAMES = "Center,Verification_Status,Validation_Status,Mutation_Status,Sequencing_Phase,Sequence_Source,Validation_Method,Score,BAM_file,Sequencer,Tumor_Sample_UUID,Match_Norm_Sample_UUID,Caller";
my $output = dirname($vcf);
# ## Create tmp and ref directory. Delete these later*
 if( ! -d "$output/tmp_$somatic/" ){
     print "$output/tmp_$somatic/ does not exist. Will create it now\n";
         mkdir("$output/tmp_$somatic", 0755) or die "Making tmp_$somatic didn't work $!";
}
if( ! -d "$output/ref_$somatic/" ){
    print "$output/ref_$somatic/ does not exist. Will create it now\n";
    mkdir("$output/ref_$somatic", 0755) or die "Making ref_$somatic didn't work $!";
}

my $ref_base = basename($REF_FASTA);

# softlink reference
symlink($REF_FASTA, "$output/ref_$somatic/$ref_base");
symlink("$REF_FASTA.fai", "$output/ref$somatic/$ref_base.fai");


## vep-forks is at 4 because that is how many CPUs we ask for
#print "\n/opt/common/CentOS_6/bin/v1/perl /opt/common/CentOS_6/vcf2maf/v1.6.1/maf2maf.pl --tmp-dir $output/tmp_$somatic --ref-fasta $output/ref_$somatic/$ref_base --vep-forks 4 --species $VEP_SPECIES --vep-path $VEP --vep-data $VEP --retain-cols $VEP_COLUMN_NAMES --input-maf $vcf\_$somatic\_maf1.txt --output-maf $vcf\_$somatic\_maf2_VEP.txt --ncbi-build $NCBI_BUILD \n\n";

#`/opt/common/CentOS_6/bin/v1/perl /opt/common/CentOS_6/vcf2maf/v1.6.1/maf2maf.pl --tmp-dir $output/tmp_$somatic --ref-fasta $output/ref_$somatic/$ref_base --vep-forks 4 --vep-path $VEP --species $VEP_SPECIES --vep-data $VEP --retain-cols $VEP_COLUMN_NAMES --input-maf $vcf\_$somatic\_maf1.txt --output-maf $vcf\_$somatic\_maf2_VEP.txt --ncbi-build $NCBI_BUILD > $vcf\_$somatic\_maf2_VEP.log 2>&1`;

print "/opt/common/CentOS_6/bin/v1/perl /opt/common/CentOS_6/vcf2maf/v1.5.4/maf2maf.pl --tmp-dir $output/tmp_$somatic --ref-fasta $output/ref_$somatic/$ref_base --vep-forks 4 --vep-path $VEP --vep-data $VEP --retain-cols $VEP_COLUMN_NAMES --input-maf $vcf\_$somatic\_maf1.txt --output-maf $vcf\_$somatic\_maf2_VEP.txt > $vcf\_$somatic\_maf2_VEP.log 2>&1";
                 
`/opt/common/CentOS_6/bin/v1/perl /opt/common/CentOS_6/vcf2maf/v1.5.4/maf2maf.pl --tmp-dir $output/tmp_$somatic --ref-fasta $output/ref_$somatic/$ref_base --vep-forks 4 --vep-path $VEP --vep-data $VEP --retain-cols $VEP_COLUMN_NAMES --input-maf $vcf\_$somatic\_maf1.txt --output-maf $vcf\_$somatic\_maf2_VEP.txt > $vcf\_$somatic\_maf2_VEP.log 2>&1`;

print "creating TCGA-formatted MAF file... \n";
#This removes any records that don't have a gene name at the front
`grep -v ^Unknown $vcf\_$somatic\_maf2_VEP.txt  > $vcf\_$somatic\_VEP_MAF.txt`;
`grep -v ^# $vcf\_$somatic\_VEP_MAF.txt | cut -f-34 > $vcf\_$somatic\_TCGA_MAF.txt`;

print "creating MAF for cbio portal submission";
print "$PYTHON/python $Bin/maf/pA_reSortCols.py -i $vcf\_$somatic\_VEP_MAF.txt -f $Bin/maf/finalCols_PORTAL.txt -o $vcf\_$somatic\_TCGA_PORTAL_MAF.txt\n\n";
`$PYTHON/python $Bin/maf/pA_reSortCols.py -i $vcf\_$somatic\_VEP_MAF.txt -f $Bin/maf/finalCols_PORTAL.txt -o $vcf\_$somatic\_TCGA_PORTAL_MAF.txt`;

print "creating clean MAF file\n";
print "$PYTHON/python $Bin/maf/pA_reSortCols.py -i $vcf\_$somatic\_VEP_MAF.txt -f $Bin/maf/finalCols.txt -o $vcf\_$somatic\_MAF4.txt\n\n";
# Create nice MAF with essential columns
`$PYTHON/python $Bin/maf/pA_reSortCols.py -i $vcf\_$somatic\_VEP_MAF.txt -f $Bin/maf/finalCols.txt -o $vcf\_$somatic\_MAF4.txt`;

print "annotating with cosmic\n";
print "$PYTHON/python $Bin/maf/maf_annotations/addCosmicAnnotation.py -i $vcf\_$somatic\_VEP_MAF.txt -o $vcf\_$somatic\_VEP_COSMIC_MAF_STANDARD.txt -f $Bin/data/CosmicMutantExport_v67_241013.tsv\n\n";
print "$PYTHON/python $Bin/maf/maf_annotations/addCosmicAnnotation.py -i $vcf\_$somatic\_VEP_MAF.txt -o $vcf\_$somatic\_VEP_COSMIC_MAF_DETAILED.txt -f $Bin/data/CosmicMutantExport_v67_241013.tsv -d\n\n";
`$PYTHON/python $Bin/maf/maf_annotations/addCosmicAnnotation.py -i $vcf\_$somatic\_VEP_MAF.txt -o $vcf\_$somatic\_VEP_COSMIC_MAF_STANDARD.txt -f $Bin/data/CosmicMutantExport_v67_241013.tsv`;
`$PYTHON/python $Bin/maf/maf_annotations/addCosmicAnnotation.py -i $vcf\_$somatic\_VEP_MAF.txt -o $vcf\_$somatic\_VEP_COSMIC_MAF_DETAILED.txt -f $Bin/data/CosmicMutantExport_v67_241013.tsv -d`;

print "annotating with mutation assessor\n";
print "$PYTHON/python $Bin/maf/maf_annotations/addMAannotation.py -i $vcf\_$somatic\_VEP_COSMIC_MAF_STANDARD.txt -o $vcf\_$somatic\_VEP_COSMIC_MA_MAF_STANDARD.txt\n";
print "$PYTHON/python $Bin/maf/maf_annotations/addMAannotation.py -i $vcf\_$somatic\_VEP_COSMIC_MAF_DETAILED.txt -o $vcf\_$somatic\_VEP_COSMIC_MA_MAF_DETAILED.txt -d\n\n";
`$PYTHON/python $Bin/maf/maf_annotations/addMAannotation.py -i $vcf\_$somatic\_VEP_COSMIC_MAF_STANDARD.txt -o $vcf\_$somatic\_VEP_COSMIC_MA_MAF_STANDARD.txt`;
`$PYTHON/python $Bin/maf/maf_annotations/addMAannotation.py -i $vcf\_$somatic\_VEP_COSMIC_MAF_DETAILED.txt -o $vcf\_$somatic\_VEP_COSMIC_MA_MAF_DETAILED.txt -d`;

if($patient && $bam_dir){
    #if($pairing){

        print "Starting maf fillout\n";
        # open patient file, get each sample name:
        # then find file with that name in the alignement directory
        # make sure it there is only 1 bam per sample
        # add that to a array
        open(PATIENT, "$patient") || die "Can't open patient file $patient $!";
        my $header = <PATIENT>;
        my @header = split(/\s+/,$header);

        my ($sID_index) = grep {$header[$_] =~ /Sample_ID/} 0..$#header;
        #print "Sample index: $sID_index\n";

        my @bamList;
 
        while(<PATIENT>) {
            chomp;
            my @patient=split(/\s+/,$_);
            #print "Sample: $patient[$sID_index] \n";
            my $bamFile = `find $bam_dir -name "Proj_*_indelRealigned_recal_$patient[$sID_index].bam"`;
            chomp($bamFile);
            #print "Bam file: $bamFile \n";
        
            push(@bamList, "--bam $patient[$sID_index]:$bamFile");
        }

        my $bam_inputs = join(" ", @bamList);

        print "$Bin/maf/fillout/GetBaseCountsMutliSample/GetBaseCountsMultiSample --fasta $REF_FASTA $bam_inputs --output $vcf\_$somatic\_TCGA_basecounts.txt --maf $vcf\_$somatic\_TCGA_PORTAL_MAF.txt --filter_improper_pair 0\n\n";
        `$Bin/maf/fillout/GetBaseCountsMutliSample/GetBaseCountsMultiSample --fasta $REF_FASTA $bam_inputs --output $vcf\_$somatic\_TCGA_basecounts.txt --maf $vcf\_$somatic\_TCGA_PORTAL_MAF.txt --filter_improper_pair 0`;

    if($pairing){
        print "$PYTHON/python $Bin/maf/fillout/dmp2portalMAF -m $vcf\_$somatic\_TCGA_PORTAL_MAF.txt -p $pairing -P $patient -c $caller -b $vcf\_$somatic\_TCGA_basecounts.txt -o $vcf\_$somatic\_TCGA_PORTAL_MAF_fillout.txt\n";    
        `$PYTHON/python $Bin/maf/fillout/dmp2portalMAF -m $vcf\_$somatic\_TCGA_MAF.txt -p $pairing -P $patient -c $caller -b $vcf\_$somatic\_TCGA_basecounts.txt -o $vcf\_$somatic\_TCGA_PORTAL_MAF_fillout.txt`;
    } else {
        print "$PYTHON/python $Bin/maf/fillout/dmp2portalMAF -m $vcf\_$somatic\_TCGA_PORTAL_MAF.txt -P $patient -c $caller -b $vcf\_$somatic\_TCGA_basecounts.txt -o $vcf\_$somatic\_TCGA_PORTAL_MAF_fillout.txt\n";
        `$PYTHON/python $Bin/maf/fillout/dmp2portalMAF -m $vcf\_$somatic\_TCGA_MAF.txt -P $patient -c $caller -b $vcf\_$somatic\_TCGA_basecounts.txt -o $vcf\_$somatic\_TCGA_PORTAL_MAF_fillout.txt`;
    }
}


if($delete_temp){
    `/bin/rm $vcf\_PAIRED.maf $vcf\_$somatic\_maf0.log $vcf\_UNPAIRED.maf $vcf\_$somatic\_UNFILTERED.txt $vcf\_$somatic\_maf2.txt $vcf\_$somatic\_maf3.txt $vcf\_$somatic\_MAF3_HUGO.log $vcf\_$somatic\_maf3.txt_ambiguous $vcf\_$somatic\_maf3.txt_hugo_modified $vcf\_$somatic\_VEP_COSMIC_MAF_STANDARD.txt $vcf\_$somatic\_VEP_COSMIC_MAF_DETAILED.txt $vcf\_$somatic\_UNFILTERED.txt $vcf\_$somatic\_rescued.txt $vcf\_$somatic\_maf0_rescue.log`;
    `/bin/rm -r $output/tmp_$somatic $output/ref_$somatic`;
}
