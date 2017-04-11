#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin); 
use File::Basename;

## INPUT: VCF file
## OUTPUT: MAF Files

my $VEP_COLUMN_NAMES = "Center,Verification_Status,Validation_Status,Mutation_Status,Sequencing_Phase,Sequence_Source,Validation_Method,Score,BAM_file,Sequencer,Tumor_Sample_UUID,Matched_Norm_Sample_UUID,Caller";

## to help with selective restarting
my $force_run;

my ($vcf, $pairing, $patient, $bam_dir, $species, $config, $caller, $normal_sample, $tumor_sample, $delete_temp, $exac_vcf);
GetOptions ('vcf=s' => \$vcf,
	    'species=s' => \$species,
	    'config=s' => \$config,
	    'caller=s' => \$caller,
            'align_dir=s' => \$bam_dir,
            'patient=s' => \$patient,
            'exac_vcf=s' => \$exac_vcf,
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

if($exac_vcf){
    if( !-e $exac_vcf){
        die "$exac_vcf DOES NOT EXIST";
    }
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

if($species !~ /b37|hg19|mm10|mouse/i){
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
my $VCFTOOLS = '';
my $MM10_FASTA = '';
my $MM10_CUSTOM_FASTA = '';
my $B37_FASTA = '';
my $B37_MM10_HYBRID_FASTA = '';
my $FIXMULTIINDEL = '';

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
    elsif($conf[0] =~ /mm10_custom_fasta/i){
        if(!-e "$conf[1]"){
            if($species =~ /mm10_custom/i){
                die "CAN'T FIND $conf[1] $!";
            }
        }
        $MM10_CUSTOM_FASTA = $conf[1];
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
    elsif($conf[0] =~/vcftools/i){
        if(!-e "$conf[1]/vcftools"){
            die "CAN'T FIND vcftools in $conf[1] $!";
        }
        $VCFTOOLS = $conf[1];
    }
    elsif($conf[0] =~ /bcftools/i){
        if(!-e "$conf[1]/bcftools"){
            die "CAN'T FIND bcftools IN $conf[1] $!";
        }
        my $path_tmp = $ENV{'PATH'};
        $ENV{'PATH'} = "$conf[1]:$path_tmp";
    }
    elsif($conf[0] =~ /bedtools/i){
        if(!-e "$conf[1]/bedtools"){
            die "CAN'T FIND bedtools IN $conf[1] $!";
        }
        my $path_tmp = $ENV{'PATH'};
        $ENV{'PATH'} = "$conf[1]:$path_tmp";
    }
    elsif($conf[0] =~/fixmultiindel/i){
        if(!-e "$conf[1]/fixMultiInDel.sh"){
            die "CAN'T FIND fixMultiInDel in $conf[1] $!";
        }
        $FIXMULTIINDEL = $conf[1];
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

if(!-e $VCFTOOLS){
    die "Cannot find vcftools. This is necessary to the script. $!";
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
} # Hybrid is to be run as human, but with the same fasta it was run with before.
elsif($species =~ /hybrid/i){
    $species = 'b37';
    $REF_FASTA = "$B37_MM10_HYBRID_FASTA";
    $NCBI_BUILD = "GRCh37";
    $VEP_SPECIES = "homo_sapiens";
}
elsif($species =~ /mouse|^mm10$/i){
    $species = 'mm10';
    $REF_FASTA = "$MM10_FASTA";
    $NCBI_BUILD = "GRCm38";
    $VEP_SPECIES = "mus_musculus";
}
elsif($species =~ /^mm10_custom$/i){
    $species = 'mm10_custom';
    $REF_FASTA = "$MM10_CUSTOM_FASTA";
    $NCBI_BUILD = "GRCm38";
    $VEP_SPECIES = "mus_musculus";
}
my $output = dirname($vcf);

my $progress = "$output/progress_$somatic";

my $numLines = `grep -c -v "^#" $vcf`;
if($numLines == 0){
    `touch $vcf\_UNPAIRED_TCGA_MAF.txt $vcf\_UNPAIRED_TCGA_PORTAL_MAF.txt $vcf\_UNPAIRED_TCGA_PORTAL_MAF_fillout.txt $vcf\_UNPAIRED_VEP_MAF.txt `;
    exit 0;
}

if (! -d "$progress"){
    mkdir("$progress", 0755) or die "Making progress directory did not work. $!";
}

###
### Converting to MAF
###

my @indv_mafs;
#Haplotype caller
if($caller =~ /unifiedgenotyper|ug|haplotypecaller|hc/i){
    ## step zero. Fix multi Indels:
    print "$FIXMULTIINDEL/fixMultiInDel.sh -a $vcf $vcf\_split.vcf";
    `$FIXMULTIINDEL/fixMultiInDel.sh -a $vcf $vcf\_split.vcf`;


    ## step one, create a lovely script to add "REF.<species>" sample into the vcf so that if there is an unpaired sample, or all of the sample are unpaired, we can still do vcf2maf.
    if($force_run || !-e "$progress/" . basename("$vcf") . "_add_REF.done"){
        $force_run = 1;
        my $outVCF = "$vcf\_REF.vcf";
        print "$PERL/perl $Bin/maf/add_ref_sample_to_vcf.pl -in_vcf $vcf\_split.vcf -out $outVCF -species $species\n";
        `$PERL/perl $Bin/maf/add_ref_sample_to_vcf.pl -in_vcf $vcf\_split.vcf -out $outVCF -species $species`;
        &checkResult($?, $progress, basename("$vcf") . "_add_REF");
    }
    
    # To use vcf2maf we first have to split the vcf to pairs. This will be output in $output/tmp_$somatic
    if( ! -d "$output/tmp_$somatic/" ){
        print "$output/tmp_$somatic/ does not exist. Will create it now\n";
        mkdir("$output/tmp_$somatic", 0755) or die "Making tmp_$somatic didn't work $!";
    }
    
    my $refNormal_id =  "REF." . $species;
    
    #Split vcf here. Then call vcf2maf for each split
    if($pairing){
        open(PAIRS, $pairing);
        while(my $line = <PAIRS>){
            chomp $line;
            my @splitLine = split(/\s+/, $line);
            my $tumor = $splitLine[0];
            my $normal = $splitLine[1];
           
            if($normal =~ /na/i || $tumor =~ /na/i){
                next;
            }
            my $split_vcf_prefix = "$output/tmp_$somatic/" . basename($vcf) . "_$tumor\_$normal";
                 
            ## split vcf into paired vcfs       
            if($force_run || !-e "$progress/" . basename("$vcf") . "_$tumor\_$normal\_vcf.done"){
                $force_run = 1;
                &splitVcf("$vcf\_REF.vcf", $split_vcf_prefix, $tumor, $normal, "$progress/",  basename("$vcf") . "_$tumor\_$normal\_vcf", $split_vcf_prefix . ".recode.vcf");    
            }
            
            ## Filtering for actual variants (since we are splitting a vcf and the vcftools puts all records in vcf output)
            if($force_run || !-e "$progress/" . basename("$vcf") . "_$tumor\_$normal\_vcf_filterNonVars.done"){
                $force_run = 1;
                print "$PYTHON/python $Bin/maf/filter_nonVars_vcf.py -v $split_vcf_prefix\.recode.vcf -t $tumor -n $normal --somatic -o $split_vcf_prefix\_vars.vcf \n";
                `$PYTHON/python $Bin/maf/filter_nonVars_vcf.py -v $split_vcf_prefix\.recode.vcf -t $tumor -n $normal --somatic -o $split_vcf_prefix\_vars.vcf`;
                &checkResult($?, $progress, basename("$vcf") . "_$tumor\_$normal\_vcf_filterNonVars")
            }
            
            ## run vcf2maf ONLY if sample's vcf is not empty (or one line)
            my $numLinesVCF = `grep -c -v "^#" $split_vcf_prefix\_vars.vcf`;
            if($numLinesVCF > 0){
                if($force_run || !-e "$progress/" . basename("$vcf") . "_$tumor\_$normal\_vcf2maf.done"){
                    $force_run = 1;
                    &run_vcf_2_maf("$split_vcf_prefix\_vars.vcf", $split_vcf_prefix . ".maf", $tumor, $normal,  basename("$vcf") . "_$tumor\_$normal\_vcf2maf");
                }
                push(@indv_mafs, $split_vcf_prefix . ".maf");
            }
        }
    }else{
        ## How do I find all the sample names ? :(
        my $ref_vcf = "$vcf\_REF.vcf";
        my $header_part = `grep ^#CHROM $ref_vcf | cut -f10- `;
        my @samples = split(/\s+/, $header_part);
        my $unpaired_norm = $refNormal_id;
        for my $x (@samples){
            if($x =~ "$unpaired_norm"){
                next;
            }
            my $tumor = $x;
            my $split_vcf_prefix = "$output/tmp_$somatic/" . basename($vcf) . "_$tumor\_$unpaired_norm";
            
            ## split vcf into "paired" vcfs
            if($force_run || !-e "$progress/" . basename("$vcf") . "_$tumor\_$unpaired_norm\_vcf.done"){
                $force_run = 1;
                &splitVcf("$vcf\_REF.vcf", $split_vcf_prefix, $tumor, $unpaired_norm, "$progress/",  basename("$vcf") . "_$tumor\_$unpaired_norm\_vcf", $split_vcf_prefix . ".recode.vcf");    
            }
            
            ## Filtering for actual variants (since we are splitting a vcf and the vcftools puts all records in vcf output)
            if($force_run || !-e "$progress/" . basename("$vcf") . "_$tumor\_$unpaired_norm\_vcf_filterNonVars.done"){
                $force_run = 1;
                print "$PYTHON/python $Bin/maf/filter_nonVars_vcf.py -v $split_vcf_prefix\.recode.vcf -t $tumor -n $unpaired_norm -o $split_vcf_prefix\_vars.vcf \n";
                `$PYTHON/python $Bin/maf/filter_nonVars_vcf.py -v $split_vcf_prefix\.recode.vcf -t $tumor -n $unpaired_norm -o $split_vcf_prefix\_vars.vcf`;
                &checkResult($?, $progress, basename("$vcf") . "_$tumor\_$unpaired_norm\_vcf_filterNonVars")
            }
            
            ## run vcf2ma ONLY if sample's vcf is not empty (or one line)
            my $numLinesVCF = `grep -c -v "^#" $split_vcf_prefix\_vars.vcf`;
            if($numLinesVCF > 0){
                if($force_run || !-e "$progress/" . basename("$vcf") . "_$tumor\_$unpaired_norm\_vcf2maf.done"){
                    $force_run = 1;
                    &run_vcf_2_maf("$split_vcf_prefix\_vars.vcf", $split_vcf_prefix . ".maf", $tumor, $unpaired_norm,  basename("$vcf") . "_$tumor\_$unpaired_norm\_vcf2maf");
                }
            
                push(@indv_mafs, $split_vcf_prefix . ".maf");
            }
        }
    }
       
}

#COMBINE MAF FILES
my $raw_maf = " $output/" . basename("$vcf") . "_RAW_MAF.txt";
if($force_run || !-e "$progress/" . basename("$vcf") . "_combine_MAF.done"){
    $force_run = 1;
    # Take array of all tumor/normal pairs, and add them together.
    my $mafs = join(" ", @indv_mafs);
    
    print 'grep -e ^# -e ^Hugo_Symbol $indv_mafs[0] > $raw_maf\n';
    `grep -e ^# -e ^Hugo_Symbol $indv_mafs[0] > $raw_maf`;
    
    print "cat $mafs | grep -v -e ^# -e ^Hugo_Symbol >> $raw_maf\n";
    `cat $mafs | grep -v -e ^# -e ^Hugo_Symbol >> $raw_maf`;

    #&checkResult($?, $progress,  basename("$vcf") . "_combine_MAF");
}

# check line numbers after vcf2maf
$numLines = `grep -c -v "^#" $raw_maf`;
#if num lines == 1 (header line)
if($numLines == 1){
    `touch $vcf\_UNPAIRED_TCGA_MAF.txt $vcf\_UNPAIRED_TCGA_PORTAL_MAF.txt $vcf\_UNPAIRED_TCGA_PORTAL_MAF_fillout.txt $vcf\_UNPAIRED_VEP_MAF.txt `;
    exit 0;
}


###
###  Finished converting to MAF.
###


##  BIC general filters (LowQual and depth/freq depending on pairing)

if($caller =~ /unifiedgenotyper|ug|haplotypecaller|hc/i){
    if($force_run || ! -e "$progress/" . basename("$vcf") . "_pA_qFiltersHC.done"){
        $force_run=1;
        my $somaticPart = "";
        if($pairing){
            $somaticPart = "-s";
        }
        
        print "$PYTHON/python $Bin/maf/pA_qFiltersHC.py -m $raw_maf $somaticPart -c $caller -o $vcf\_BIC.maf";
        `$PYTHON/python $Bin/maf/pA_qFiltersHC.py -m $raw_maf $somaticPart -c $caller -o $vcf\_BIC.maf`;
        &checkResult($?, $progress, basename("$vcf") . "_pA_qFiltersHC");
    }
}


#########
##
#
#
#   BLAH!
#
#
##
#########

# Check line nums again
$numLines = `grep -v "^#" $vcf\_BIC.maf | grep -c -v "^Hugo_Symbol"`;
if($numLines == 0){
    unlink("$vcf\_BIC.maf");
    `touch $vcf\_UNPAIRED_TCGA_MAF.txt $vcf\_UNPAIRED_TCGA_PORTAL_MAF.txt $vcf\_UNPAIRED_TCGA_PORTAL_MAF_fillout.txt $vcf\_UNPAIRED_VEP_MAF.txt `;
    exit 0;
}
        

if($force_run || !-e "$progress/" . basename("$vcf") . "_TCGA_MAF.done"){
    $force_run = 1;
    print "creating TCGA-formatted MAF file... \n";
    #This removes any records that don't have a gene name at the front
    #`grep -v ^Unknown $vcf\_$somatic\_maf1.VEP  > $vcf\_$somatic\_VEP_MAF.txt`;
    `cut -f-34 $vcf\_BIC.maf > $vcf\_$somatic\_TCGA_MAF.txt`;
    &checkResult($?, $progress, basename("$vcf") . "_TCGA_MAF");
}


if($force_run || !-e "$progress/" . basename("$vcf") . "_pA_resortCols.done"){
    $force_run = 1;
    print "creating MAF for cbio portal submission";
    print "$PYTHON/python $Bin/maf/pA_reSortCols.py -i $vcf\_BIC.maf -f $Bin/maf/finalCols_PORTAL.txt -o $vcf\_$somatic\_TCGA_PORTAL_MAF.txt\n\n";
    `$PYTHON/python $Bin/maf/pA_reSortCols.py -i $vcf\_BIC.maf  -f $Bin/maf/finalCols_PORTAL.txt -o $vcf\_$somatic\_TCGA_PORTAL_MAF.txt`;
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
            my $bamFile = `find -L $bam_dir -name "*_indelRealigned_recal_$patient[$sID_index].bam"`;
            chomp($bamFile);
            #print "Bam file: $bamFile \n";
        
            push(@bamList, "--bam $patient[$sID_index]:$bamFile");
        }

        my $bam_inputs = join(" ", @bamList);

        print "$Bin/maf/fillout/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta $REF_FASTA $bam_inputs --output $vcf\_$somatic\_TCGA_basecounts.txt --maf $vcf\_$somatic\_TCGA_PORTAL_MAF.txt --filter_improper_pair 0\n\n";
        `$Bin/maf/fillout/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta $REF_FASTA $bam_inputs --output $vcf\_$somatic\_TCGA_basecounts.txt --maf $vcf\_$somatic\_TCGA_PORTAL_MAF.txt --filter_improper_pair 0 > $vcf\_$somatic\_basecounts.log 2>&1`;
        &checkResult($?, $progress, basename("$vcf") . "_getBaseCounts", "$vcf\_$somatic\_TCGA_basecounts.txt");
    } 

    my $callerString = lc $caller;

    if( !-e "$progress/" . basename("$vcf") . "_dmp2portal.done"){
        if($pairing){
            print "$PYTHON/python $Bin/maf/fillout/dmp2portalMAF -s $species -m $vcf\_$somatic\_TCGA_PORTAL_MAF.txt -p $pairing -P $patient -c $callerString -b $vcf\_$somatic\_TCGA_basecounts.txt -o $vcf\_$somatic\_TCGA_PORTAL_MAF_fillout.txt\n";    
            `$PYTHON/python $Bin/maf/fillout/dmp2portalMAF -s $species -m $vcf\_$somatic\_TCGA_PORTAL_MAF.txt -p $pairing -P $patient -c $callerString -b $vcf\_$somatic\_TCGA_basecounts.txt -o $vcf\_$somatic\_TCGA_PORTAL_MAF_fillout.txt`;
            &checkResult($?, $progress, basename("$vcf") . "_dmp2portal");
        } else {
            print "$PYTHON/python $Bin/maf/fillout/dmp2portalMAF -s $species -m $vcf\_$somatic\_TCGA_PORTAL_MAF.txt -P $patient -c $callerString -b $vcf\_$somatic\_TCGA_basecounts.txt -o $vcf\_$somatic\_TCGA_PORTAL_MAF_fillout.txt\n";
            `$PYTHON/python $Bin/maf/fillout/dmp2portalMAF -s $species -m $vcf\_$somatic\_TCGA_PORTAL_MAF.txt -P $patient -c $callerString -b $vcf\_$somatic\_TCGA_basecounts.txt -o $vcf\_$somatic\_TCGA_PORTAL_MAF_fillout.txt`;
            &checkResult($?, $progress, basename("$vcf") . "_dmp2portal");
        }
    }
}


if($species =~ /hg19|human|b37/){
    if($force_run || !-e "$progress/" . basename("$vcf") . "_bedtools_anno.done"){
        $force_run = 1;
        print "$PERL/perl $Bin/maf/bedtools_annotations.pl --in_maf $vcf\_BIC.maf --species $species --output $output --config $config --fastq --target $Bin/targets/IMPACT410_$species/IMPACT410_$species\_targets_plus5bp.bed --targetname IMPACT_410 --somatic $somatic\n";
        `$PERL/perl $Bin/maf/bedtools_annotations.pl --in_maf $vcf\_BIC.maf --species $species --output $output --config $config --fastq --target $Bin/targets/IMPACT410_$species/IMPACT410_$species\_targets_plus5bp.bed --targetname IMPACT_410 --somatic $somatic`;
        &checkResult($?, $progress, basename("$vcf") . "_bedtools_anno");
    }

    # vcf2maf includes exac annotation
    #if($force_run || !-e "$progress/" . basename("$vcf") . "_exac_anno.done"){
    #    $force_run = 1;
    #    print "perl $Bin/maf/exac_annotate.pl --in_maf $vcf\_BIC.maf --species $species --output $output --config $config --somatic $somatic --data $Bin/data\n";
    #    `$PERL/perl $Bin/maf/exac_annotate.pl --in_maf $vcf\_BIC.maf --species $species --output $output --config $config --somatic $somatic --data $Bin/data`;
    #    &checkResult($?, $progress, basename("$vcf") . "_exac_anno");
    #}

    if($force_run || !-e "$progress/" . basename("$vcf") . "_mergeExtraCols.done"){
        $force_run = 1;
        print "$PYTHON/python $Bin/maf/mergeExtraCols.py --triNuc $output/TriNuc.txt --targeted $output/maf_targets.IMPACT410 --maf $vcf\_BIC.maf\n";
       `$PYTHON/python $Bin/maf/mergeExtraCols.py --triNuc $output/TriNuc.txt --targeted $output/maf_targets.IMPACT_410 --maf $vcf\_BIC.maf > $vcf\_$somatic\_VEP_MAF.txt`;
        &checkResult($?, $progress, basename("$vcf") . "_mergeExtraCols");
    }
}else{ ## MOUSE
    if($force_run || !-e "$progress/" . basename("$vcf") . "_bedtools_anno.done"){
        $force_run = 1;
        print "$PERL/perl $Bin/maf/bedtools_annotations.pl --in_maf $vcf\_BIC.maf --species $species --output $output --config $config --fastq --somatic $somatic\n";
        `$PERL/perl $Bin/maf/bedtools_annotations.pl --in_maf $vcf\_BIC.maf --species $species --output $output --config $config --fastq --somatic $somatic`;
        &checkResult($?, $progress, basename("$vcf") . "_bedtools_anno");
    }

    if($force_run || !-e "$progress/" . basename("$vcf") . "_mergeExtraCols.done"){
        $force_run = 1;
        `/bin/touch $output/blank`;
        print "$PYTHON/python $Bin/maf/mergeExtraCols.py --triNuc $output/TriNuc.txt --maf $vcf\_BIC.maf > $vcf\_$somatic\_VEP_MAF.txt\n";
        `$PYTHON/python $Bin/maf/mergeExtraCols.py --triNuc $output/TriNuc.txt --maf $vcf\_BIC.maf > $vcf\_$somatic\_VEP_MAF.txt`;
        &checkResult($?, $progress, basename("$vcf") . "_mergeExtraCols", "$vcf\_$somatic\_VEP_MAF.txt");
    }
}


#####
#
#
#

#
#
#
#
#
#####



sub run_vcf_2_maf{
    my ($in_vcf, $out_maf, $tumor, $normal, $doneFile) = @_;
    
    if( ! -d "$output/ref_$somatic/" ){
        print "$output/ref_$somatic/ does not exist. Will create it now\n";
        mkdir("$output/ref_$somatic", 0755) or die "Making ref_$somatic didn't work $!";
    }

    my $ref_base = basename($REF_FASTA);
    
    # softlink reference
    symlink($REF_FASTA, "$output/ref_$somatic/$ref_base");
    symlink("$REF_FASTA.fai", "$output/ref_$somatic/$ref_base.fai");
    
    my $addOptions = '';
    if($exac_vcf){
        $addOptions = "--filter-vcf $exac_vcf";
    }
    
    print "\n#######\n#######\nStarting VEP. \n";
    print "$PERL/perl $VCF2MAF/vcf2maf.pl --input-vcf $in_vcf --species $VEP_SPECIES --output-maf $out_maf --ref-fasta $output/ref_$somatic/$ref_base --tmp-dir $output/tmp_$somatic/ --ncbi $NCBI_BUILD --vep-forks 4 --vep-path $VEP --vep-data $VEP $addOptions --tumor-id $tumor --normal-id $normal \n";
    `$PERL/perl $VCF2MAF/vcf2maf.pl --input-vcf $in_vcf --species $VEP_SPECIES --output-maf $out_maf --ref-fasta $output/ref_$somatic/$ref_base --tmp-dir $output/tmp_$somatic/ --ncbi $NCBI_BUILD --vep-forks 4 --vep-path $VEP --vep-data $VEP $addOptions --tumor-id $tumor --normal-id $normal`;

    &checkResult($?, $progress, $doneFile, $out_maf);
}

sub splitVcf{
    my ($inVCF, $outPrefix, $tumor, $normal, $progress, $doneFile, $checkFile) = @_;
    #Now I have to use vcftools to separate the paired stuff
    print "$VCFTOOLS/vcftools --indv $tumor --indv $normal --vcf $inVCF --out $outPrefix --recode\n";
    `$VCFTOOLS/vcftools --indv $tumor --indv $normal --vcf $inVCF --out $outPrefix --recode`;
    &checkResult($?, $progress, $doneFile, $checkFile);
}

sub checkResult{
    my ($status, $progress, $filebase, $out_check) = @_;

    if($out_check){
        if($status == 0 && -e "$out_check" && ! -z "$out_check" ) 
        { 
            `/bin/touch $progress/$filebase.done`;
        }
    } elsif ($status == 0){
        `/bin/touch $progress/$filebase.done`;
    } else{
        exit "\nThere was an error with script, or the output file was not created $filebase!\n";
    }
    

}


