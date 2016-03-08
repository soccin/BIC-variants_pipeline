#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin); 
use File::Basename;

### INPUT: vcf file and list of normal/tumor sample pairing information
### OUTPUT: 2 maf files; 

## CONSTANT FOR VEP
my $VEP_COLUMN_NAMES = "Center,Verification_Status,Validation_Status,Mutation_Status,Sequencing_Phase,Sequence_Source,Validation_Method,Score,BAM_file,Sequencer,Tumor_Sample_UUID,Matched_Norm_Sample_UUID,Caller";

## to help with selective restarting
my $force_run;

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

if($species !~ /b37|hg19/i){ ##    UNCOMMENT LATER   |mm10|mouse/i){
    print "THIS WILL ONLY PRINT OUT A TCGA MAF, NO ANNOATION\n";
}

if($caller !~ /unifiedgenotyper|ug|haplotypecaller|hc|mutect|varscan|somaticsniper/i){
    die "Only support for unifiedgenotyper(ug), haplotypecaller(hc), mutect, varscan, and somaticsniper";
}

my $REF_FASTA = '';
my $HG19_FASTA = '';
my $PYTHON = '';
my $PERL = '';
my $VCF2MAF = '';
my $VEP = '';
my $MM10_FASTA = '';
my $B37_FASTA = '';
my $B37_MM10_HYBRID_FASTA = '';

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
    elsif($conf[0] =~ /b37_fasta/i){
        if(!-e "$conf[1]"){
            if($species =~ /human|^b37$/i){
                die "CAN'T FIND $conf[1] $!";
            }
        }
        $B37_FASTA = $conf[1];
    }
    elsif($conf[0] =~ /b37_mm10_hybrid_fasta/i){
        if(!-e "$conf[1]"){
            if($species =~ /hybrid|b37_mm10/i){
                die "CAN'T FIND $conf[1] $!";
            }
        }
        $B37_MM10_HYBRID_FASTA = $conf[1];
    }
    elsif($conf[0] =~ /^perl/i){
        if(!-e "$conf[1]/perl"){
            die "CAN'T FIND perl IN $conf[1] $!";
        }
        $PERL = $conf[1];
    }
    elsif($conf[0] =~ /^vep/i){
        if(!-e "$conf[1]/variant_effect_predictor.pl"){
            die "CAN'T FIND VEP IN $conf[1] $!";
        }
        $VEP = $conf[1];
    }
    elsif($conf[0] =~/vcf2maf/i){
        if(!-e "$conf[1]/maf2maf.pl"){
            die "CAN'T FIND maf2maf.pl in $conf[1] $!";
        }
        $VCF2MAF = $conf[1];
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

if(!-e $VEP ){
    die "VEP path from config file does not exist";
}
if(!-e $VCF2MAF ) {
    die "Cannot find vcf2maf. Either add or correct config file. $!";
}

my $NCBI_BUILD = '';
my $VEP_SPECIES = '';
if($species =~ /hg19/i){
    $species = 'hg19';
    $REF_FASTA = "$HG19_FASTA";
    $NCBI_BUILD = "GRCh37";
    $VEP_SPECIES = "homo_sapiens";
}
elsif($species =~ /human|^b37/i){
    $species = 'b37';
    $REF_FASTA = "$B37_FASTA";
    $NCBI_BUILD = "GRCh37";
    $VEP_SPECIES = "homo_sapiens";
}
elsif($species =~ /mouse|^mm10$/i){
    $species = 'mm10';
    $REF_FASTA = "$MM10_FASTA";
    $NCBI_BUILD = "GRCm38";
    $VEP_SPECIES = "mus_musculus";
}

my $output = dirname($vcf);

my $progress = "$output/progress_$somatic";
if (! -d "$progress"){
    mkdir("$progress", 0755) or die "Making progress directory did not work. $!";
}

print "converting to MAF and filtering snp calls for coverage/qual\n";
if($caller =~ /unifiedgenotyper|ug|haplotypecaller|hc/i){
    if($pairing){
        if($force_run || !-e "$progress/" . basename("$vcf") . "_vcf2maf0.done"){ 
            $force_run = 1;
            print "$PYTHON/python $Bin/maf/vcf2maf0.py -i $vcf -c $caller -o $vcf\_PAIRED.maf -p $pairing\n";
	    `$PYTHON/python $Bin/maf/vcf2maf0.py -i $vcf -c $caller -o $vcf\_PAIRED.maf -p $pairing`;
            &checkResult($?, $progress, basename("$vcf") . "_vcf2maf0"); 
        }
        if($force_run || !-e "$progress/" . basename("$vcf") . "_pA_qsomHC.done"){
            $force_run = 1;
            print "cat $vcf\_PAIRED.maf | $PYTHON/python $Bin/maf/pA_qSomHC.py | /usr/bin/tee $vcf\_$somatic\_maf0.txt\n\n";
	    `cat $vcf\_PAIRED.maf | $PYTHON/python $Bin/maf/pA_qSomHC.py | /usr/bin/tee $vcf\_$somatic\_maf0.txt > $vcf\_$somatic\_maf0.log 2>&1`;
            &checkResult($?, $progress, basename("$vcf") . "_pA_qsomHC");
        }
    }
    else{
        if($force_run || !-e "$progress/" . basename("$vcf") . "_vcf2maf0.done"){ 
            $force_run = 1;
            print "$PYTHON/python $Bin/maf/vcf2maf0.py -i $vcf -c $caller -o $vcf\_UNPAIRED.maf\n";
	    `$PYTHON/python $Bin/maf/vcf2maf0.py -i $vcf -c $caller -o $vcf\_UNPAIRED.maf`;
            &checkResult($?, $progress, basename("$vcf") . "_vcf2maf0");
        }
        if($force_run || !-e "$progress/" . basename("$vcf") . "_pA_qsomHC.done"){
            $force_run = 1;
            print "cat $vcf\_UNPAIRED.maf | $PYTHON/python $Bin/maf/pA_Cov+noLow.py | /usr/bin/tee $vcf\_$somatic\_maf0.txt\n\n";
            `cat $vcf\_UNPAIRED.maf | $PYTHON/python $Bin/maf/pA_Cov+noLow.py | /usr/bin/tee $vcf\_$somatic\_maf0.txt > $vcf\_$somatic\_maf0.log 2>&1`;
            &checkResult($?, $progress, basename("$vcf") . "_pA_qsomHC");
        }
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

    if($force_run || !-e "$progress/" . basename("$vcf") . "_vcf2maf0.done"){
        $force_run = 1;
        print "$PYTHON/python $Bin/maf/vcf2maf0.py -i $vcf $af -c $caller -o $vcf\_$somatic\_maf0.txt -p $pairing $n_sample $t_sample\n\n";
        `$PYTHON/python $Bin/maf/vcf2maf0.py -i $vcf $af -c $caller -o $vcf\_$somatic\_maf0.txt -p $pairing $n_sample $t_sample`;
        &checkResult($?, $progress, basename("$vcf") . "_vcf2maf0");
    }

    if($caller =~ /mutect|mu/i){
        if($force_run || !-e "$progress/" . basename("$vcf") . "_rescue.done"){
            $force_run = 1;
	    `/bin/mv $vcf\_$somatic\_maf0.txt $vcf\_$somatic\_UNFILTERED.txt`;
            print "$PYTHON/python $Bin/rescue/DMP_rescue.py <$vcf\_$somatic\_UNFILTERED.txt> $vcf\_$somatic\_rescued.txt 2> $vcf\_$somatic\_maf0_rescue.log";
	    `$PYTHON/python $Bin/rescue/DMP_rescue.py <$vcf\_$somatic\_UNFILTERED.txt> $vcf\_$somatic\_rescued.txt 2> $vcf\_$somatic\_maf0_rescue.log`;
	    &checkResult($?, $progress, basename("$vcf") . "_rescue");
    
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
}

if($force_run || !-e "$progress/" . basename("$vcf") . "_oldmaf2tcgamaf.done"){
    $force_run = 1;
    print "converting to new MAF format using species $species ... \n";
    print "$PYTHON/python $Bin/maf/oldMAF2tcgaMAF.py $species \<$vcf\_$somatic\_maf0.txt\> $vcf\_$somatic\_maf1.txt\n\n";
    ### Convert old (NDS) MAF to official TCGA MAF
    ### NOTE; DON'T FORGET <> AROUND INPUT FILE (maf0.txt)
    `$PYTHON/python $Bin/maf/oldMAF2tcgaMAF.py $species $vcf\_$somatic\_maf0.txt $vcf\_$somatic\_maf1.txt`;
    &checkResult($?, $progress, basename("$vcf") . "_oldmaf2tcgamaf");
}
if($species !~ /hg19|b37|mm10|mouse|human/i) { ###   |mm10|mouse/i) { uncomment later!
    `cut -f-34 $vcf\_$somatic\_maf1.txt > $vcf\_$somatic\_TCGA_MAF.txt`;
    print "End of species ambiguous portion of the script.\n";

    if($delete_temp){
        &cleanUp;
    }
    exit 0;
}

# these are names needed for the "retain-cols" option in VEP
#/ my $VEP_COLUMN_NAMES = "Center,Verification_Status,Validation_Status,Mutation_Status,Sequencing_Phase,Sequence_Source,Validation_Method,Score,BAM_file,Sequencer,Tumor_Sample_UUID,Match_Norm_Sample_UUID,Caller";
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
symlink("$REF_FASTA.fai", "$output/ref_$somatic/$ref_base.fai");

if($force_run || !-e "$progress/" . basename("$vcf") . "_VEP.done"){
    $force_run = 1;
    print "\n#######\n#######\nStarting VEP. \n";
    print "\n$PERL/perl $VCF2MAF/maf2maf.pl --tmp-dir $output/tmp_$somatic --ref-fasta $output/ref_$somatic/$ref_base --ncbi-build $NCBI_BUILD --species $VEP_SPECIES --vep-forks 4 --vep-path $VEP --vep-data $VEP --retain-cols $VEP_COLUMN_NAMES --input-maf $vcf\_$somatic\_maf1.txt --output-maf $vcf\_$somatic\_maf1.VEP\n\n";

    `$PERL/perl $VCF2MAF/maf2maf.pl --tmp-dir $output/tmp_$somatic --ref-fasta $output/ref_$somatic/$ref_base --ncbi-build $NCBI_BUILD --species $VEP_SPECIES --vep-forks 4 --vep-path $VEP --vep-data $VEP --retain-cols $VEP_COLUMN_NAMES --input-maf $vcf\_$somatic\_maf1.txt --output-maf $vcf\_$somatic\_maf1.VEP > $vcf\_$somatic\_maf1.log 2>&1`;
     &checkResult($?, $progress, basename("$vcf") . "_VEP", "$vcf\_$somatic\_maf1.VEP");
}

if($force_run || !-e "$progress/" . basename("$vcf") . "_TCGA_MAF.done"){
    $force_run = 1;
    print "creating TCGA-formatted MAF file... \n";
    #This removes any records that don't have a gene name at the front
    #`grep -v ^Unknown $vcf\_$somatic\_maf1.VEP  > $vcf\_$somatic\_VEP_MAF.txt`;
    `cut -f-34 $vcf\_$somatic\_maf1.VEP > $vcf\_$somatic\_TCGA_MAF.txt`;
    &checkResult($?, $progress, basename("$vcf") . "_TCGA_MAF");
}


if($force_run || !-e "$progress/" . basename("$vcf") . "_pA_resortCols.done"){
    $force_run = 1;
    print "creating MAF for cbio portal submission";
    print "$PYTHON/python $Bin/maf/pA_reSortCols.py -i $vcf\_$somatic\_maf1.VEP -f $Bin/maf/finalCols_PORTAL.txt -o $vcf\_$somatic\_TCGA_PORTAL_MAF.txt\n\n";
    `$PYTHON/python $Bin/maf/pA_reSortCols.py -i $vcf\_$somatic\_maf1.VEP -f $Bin/maf/finalCols_PORTAL.txt -o $vcf\_$somatic\_TCGA_PORTAL_MAF.txt`;
    &checkResult($?, $progress, basename("$vcf") . "_pA_resortCols");
}

if($patient && $bam_dir){
    if($force_run || !-e "$progress/" . basename("$vcf") . "_getBaseCounts.done"){
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

        print "$Bin/maf/fillout/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta $REF_FASTA $bam_inputs --output $vcf\_$somatic\_TCGA_basecounts.txt --maf $vcf\_$somatic\_TCGA_PORTAL_MAF.txt --filter_improper_pair 0\n\n";
        `$Bin/maf/fillout/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta $REF_FASTA $bam_inputs --output $vcf\_$somatic\_TCGA_basecounts.txt --maf $vcf\_$somatic\_TCGA_PORTAL_MAF.txt --filter_improper_pair 0 > $vcf\_$somatic\_basecounts.log 2>&1`;
        &checkResult($?, $progress, basename("$vcf") . "_getBaseCounts", "$vcf\_$somatic\_TCGA_basecounts.txt");

        } 

    if( !-e "$progress/" . basename("$vcf") . "_dmp2portal.done"){
        if($pairing){
            print "$PYTHON/python $Bin/maf/fillout/dmp2portalMAF -s $species -m $vcf\_$somatic\_TCGA_PORTAL_MAF.txt -p $pairing -P $patient -c $caller -b $vcf\_$somatic\_TCGA_basecounts.txt -o $vcf\_$somatic\_TCGA_PORTAL_MAF_fillout.txt\n";    
            `$PYTHON/python $Bin/maf/fillout/dmp2portalMAF -s $species -m $vcf\_$somatic\_TCGA_PORTAL_MAF.txt -p $pairing -P $patient -c $caller -b $vcf\_$somatic\_TCGA_basecounts.txt -o $vcf\_$somatic\_TCGA_PORTAL_MAF_fillout.txt`;
            &checkResult($?, $progress, basename("$vcf") . "_dmp2portal");
        } else {
            print "$PYTHON/python $Bin/maf/fillout/dmp2portalMAF -s $species -m $vcf\_$somatic\_TCGA_PORTAL_MAF.txt -P $patient -c $caller -b $vcf\_$somatic\_TCGA_basecounts.txt -o $vcf\_$somatic\_TCGA_PORTAL_MAF_fillout.txt\n";
            `$PYTHON/python $Bin/maf/fillout/dmp2portalMAF -s $species -m $vcf\_$somatic\_TCGA_PORTAL_MAF.txt -P $patient -c $caller -b $vcf\_$somatic\_TCGA_basecounts.txt -o $vcf\_$somatic\_TCGA_PORTAL_MAF_fillout.txt`;
            &checkResult($?, $progress, basename("$vcf") . "_dmp2portal");
        }
    }
}


if($species =~ /hg19|human|b37/){
    if($force_run || !-e "$progress/" . basename("$vcf") . "_bedtools_anno.done"){
        $force_run = 1;
        print "$PERL/perl $Bin/maf/bedtools_annotations.pl --in_maf $vcf\_$somatic\_maf1.VEP --species $species --output $output --config $config --fastq --target $Bin/targets/IMPACT410_$species/IMPACT410_$species\_targets_plus5bp.bed --targetname impact410 --somatic $somatic\n";
        `$PERL/perl $Bin/maf/bedtools_annotations.pl --in_maf $vcf\_$somatic\_maf1.VEP --species $species --output $output --config $config --fastq --target $Bin/targets/IMPACT410_$species/IMPACT410_$species\_targets_plus5bp.bed --targetname impact410 --somatic $somatic`;
        &checkResult($?, $progress, basename("$vcf") . "_bedtools_anno");
    }

    if($force_run || !-e "$progress/" . basename("$vcf") . "_exac_anno.done"){
        $force_run = 1;
        print "perl $Bin/maf/exact_annotate.pl --in_maf $vcf\_$somatic\_maf1.VEP --species $species --output $output --config $config --somatic $somatic --data $Bin/data\n";
        `$PERL/perl $Bin/maf/exact_annotate.pl --in_maf $vcf\_$somatic\_maf1.VEP --species $species --output $output --config $config --somatic $somatic --data $Bin/data`;
        &checkResult($?, $progress, basename("$vcf") . "_exac_anno");
    }

    if($force_run || !-e "$progress/" . basename("$vcf") . "_mergeExtraCols.done"){
        $force_run = 1;
        print "$PYTHON/python $Bin/maf/mergeExtraCols.py $output/triNucleotide.seq $output/maf_targets.impact410 $output/exact.vcf $vcf\_$somatic\_maf1.VEP\n";
       `$PYTHON/python $Bin/maf/mergeExtraCols.py $output/triNucleotide.seq $output/maf_targets.impact410 $output/exact.vcf $vcf\_$somatic\_maf1.VEP > $vcf\_$somatic\_VEP_MAF.txt`;
        &checkResult($?, $progress, basename("$vcf") . "_mergeExtraCols");
    }
}else{ ## MOUSE
    if($force_run || !-e "$progress/" . basename("$vcf") . "_bedtools_anno.done"){
        $force_run = 1;
        print "$PERL/perl $Bin/maf/bedtools_annotations.pl --in_maf $vcf\_$somatic\_maf1.VEP --species $species --output $output --config $config --fastq --somatic $somatic\n";
        `$PERL/perl $Bin/maf/bedtools_annotations.pl --in_maf $vcf\_$somatic\_maf1.VEP --species $species --output $output --config $config --fastq --somatic $somatic`;
        &checkResult($?, $progress, basename("$vcf") . "_bedtools_anno");
    }

    if($force_run || !-e "$progress/" . basename("$vcf") . "_mergeExtraCols.done"){
        $force_run = 1;
        `/bin/touch $output/blank`;
        print "$PYTHON/python $Bin/maf/mergeExtraCols.py $output/triNucleotide.seq $output/blank $output/blank $vcf\_$somatic\_maf1.VEP > $vcf\_$somatic\_VEP_MAF.txt\n";
        `$PYTHON/python $Bin/maf/mergeExtraCols.py $output/triNucleotide.seq $output/blank $output/blank $vcf\_$somatic\_maf1.VEP > $vcf\_$somatic\_VEP_MAF.txt`;
        &checkResult($?, $progress, basename("$vcf") . "_mergeExtraCols", "$vcf\_$somatic\_VEP_MAF.txt");
    }
}

sub checkResult{
    my ($status, $progress, $filebase, $out_check) = @_;

    if($out_check){
        if($status == 0 && -e "$out_check" ) 
        { 
            `/bin/touch $progress/$filebase.done`;
        }
    } elsif ($status == 0){
        `/bin/touch $progress/$filebase.done`;
    } else{
        exit "\nThere was an error with script, or the output file was not created $filebase!\n";
    }
    

}

if($delete_temp){
#    `/bin/rm $vcf\_PAIRED.maf $vcf\_$somatic\_maf0.log $vcf\_$somatic\_maf0.txt $vcf\_$somatic\_maf1.txt $vcf\_$somatic\_maf1.log $vcf\_$somatic\_maf1.VEP $vcf\_UNPAIRED.maf $vcf\_$somatic\_UNFILTERED.txt $vcf\_$somatic\_UNFILTERED.txt $vcf\_$somatic\_rescued.txt $vcf\_$somatic\_maf0_rescue.log $output/exact.vcf $output/*.seq $output/maf_targets.* $output/blank $output/*basecounts.txt`;
#    `/bin/rm -rf $output/tmp_$somatic $output/ref_$somatic $output/xtra_$somatic $output/bed_$somatic $output/progress_$somatic`;
    &cleanUp;
}

sub cleanUp{
    `/bin/rm $vcf\_PAIRED.maf $vcf\_$somatic\_maf0.log $vcf\_$somatic\_maf0.txt $vcf\_$somatic\_maf1.txt $vcf\_$somatic\_maf1.log $vcf\_$somatic\_maf1.VEP $vcf\_UNPAIRED.maf $vcf\_$somatic\_UNFILTERED.txt $vcf\_$somatic\_UNFILTERED.txt $vcf\_$somatic\_rescued.txt $vcf\_$somatic\_maf0_rescue.log $output/exact.vcf $output/*.seq $output/maf_targets.* $output/blank $output/*basecounts.txt`;
    `/bin/rm -rf $output/tmp_$somatic $output/ref_$somatic $output/xtra_$somatic $output/bed_$somatic $output/progress_$somatic`;
}
