#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin);
use File::Path qw(make_path);
use File::Basename;

my ($svnRev, $pairing, $hc_vcf, $help, $mutect_dir, $species, $bam_dir, $patient, $pre, $output, $delete_temp, $config);

# Default output = results/variants/FinalReport
# get current directory
$output = "results/variants/FinalReport";
my $curDir = `pwd`;
chomp $curDir;
my $commandline = join " ", "\n", $0, @ARGV;
my @maf_header=();
my $force_run;

# Print command line.
print "$commandline\n\n";

GetOptions ('pair=s' => \$pairing,
            'hc_vcf=s' => \$hc_vcf,
            'mutect_dir=s' => \$mutect_dir,
            'align_dir=s' => \$bam_dir,
            'patient=s' => \$patient,
            'pre=s' => \$pre,
            'output=s' => \$output,
            'delete_temp' => \$delete_temp,
	    'config=s' => \$config,
            'species=s' => \$species,
            'svnRev=s' => \$svnRev,
            'help|h' => \$help ) or exit(1);

# if the basic necessities are missing, exit
if(!$hc_vcf || !$species || !$mutect_dir || !$pre || !$config || $help){
    print <<HELP;
    
    USAGE: HaploTect_merge.pl -pair Proj_XXXX_sample_pairing.txt -hc_vcf Proj_XXXX_HaplotypeCaller.vcf -mutect_dir /ifs/work/mutect_results -pre Proj_XXXX -output /ifs/work/output 
    
    -pair\tFile listing tumor/normal pairing of samples for mutect/maf conversion; if not specified, considered unpaired
    -hc_vcf\tHaplotype caller's vcf output file
    -mutect_dir\tDirectory were MuTect output is found
    -config\tConfiguration file
    -species\tSpecies
    -pre\tProject ID Ex: Proj_4500
    -patient\tPatient file (for optional fillout)
    -align_dir\tDirectory where the aligned bams are found (for optional fillout)
    -output\toutput directory where the merge intermediate and final files will be (default : results/variants/FinalReport)
    -delete_temp\tDelete the intermediate files.
    
HELP
    exit;
}

## Check species
if($species !~ /hg19|b37|mm10/i){
    print "Haplotect only works 100% with hg19/b37/mm10 for right now, otherwise we will print out the MAF before annotation and exit\n\n";
}

## Checking HC VCF
if( ! -e $hc_vcf || -z $hc_vcf){
    die "The haplotype vcf file either does not exist or is empty. ($hc_vcf)\n";
}

if($hc_vcf !~ /^\//){
    $hc_vcf = "$curDir/$hc_vcf";
}

## Checking Mutect Directory
if( ! -d $mutect_dir ){
    die "The mutect directory either does not exist or is not a directory. ($mutect_dir)\n";
}

if($mutect_dir !~ /^\//){
    $mutect_dir = "$curDir/$mutect_dir";
}

if (! -d $output ){
    ###die "The output directory either does not exist or is not a directory. ($output)\n";
    `/bin/mkdir -m 775 -p $output`;
}

if($output !~ /^\//){
    $output = "$curDir/$output";
}

my @bamList;
if($patient) {
    if(!-e $patient){
        die "$patient DOES NOT EXIST";
    }
    if(!$bam_dir){
        die "If patient file is given, you must supply alignment directory for fillout.";
    } else {
        if(!-e $bam_dir){
            die "$bam_dir DOES NOT EXIST";
        }
        open(PATIENT, "$patient") || die "Can't open patient file $patient $!";
        my $header = <PATIENT>;
        my @header = split(/\s+/,$header);

        my ($sID_index) = grep {$header[$_] =~ /Sample_ID/} 0..$#header;
        while(<PATIENT>) {
            chomp;
            my @patient=split(/\s+/,$_);
            my $bamFile = `find -L $bam_dir -name "Proj_*_indelRealigned_recal_$patient[$sID_index].bam"`;
            chomp($bamFile);

            push(@bamList, "--bam $patient[$sID_index]:$bamFile");
        }
    }

}

if($bam_dir){
    if(!-e $bam_dir){
        die "$bam_dir DOES NOT EXIST";
    }
}


## Grab mutect vcf files from mutect dir
opendir(DIR, $mutect_dir);
my @mutect_vcfs = grep(/\.vcf$/,readdir(DIR));
closedir(DIR);

## check to make sure these are files and aren't empty(not like, directories)? (is this really necessary)
for my $vcf (@mutect_vcfs){
    if( ! -f "$mutect_dir/$vcf" || -z "$mutect_dir/$vcf"){
        die "The mutect vcf file is not a regular file or is empty. ($vcf)\n";
    }
}

#check pairing file
if(! -e $pairing){
    die "Pairing file doesn't exist\n";
}
## Checks the pairing file to make sure mutect tumor/normal vcf files are present
open(my $pair, "<", $pairing);
chomp(my @pair_lines = <$pair>);
close $pair;

## To make sure there are the same number of mutect files as there are pairs
my $pairedCount = 0;
foreach my $line (@pair_lines) {
    my($normal,$tumor) = split(/\s+/, $line);
    next if ($tumor  =~ /^NA$/i || $normal  =~ /^NA$/i );
    $pairedCount ++;
    my @vcfIndex = indices("$normal\_$tumor\_mutect", @mutect_vcfs);

    if(scalar @vcfIndex != 1){
        if(scalar @vcfIndex  < 1){
            die "No vcf file matches with these tumor/normal pair: $tumor $normal\n";
        }
        else {
            die "More than one vcf file matches with these tumor/normal pair: $tumor $normal?\n";
        }
    }
}

my $l = scalar @mutect_vcfs;
if($pairedCount != $l) {
    die "The number of pairs ($pairedCount) is not equal to the number of vcf files ($l)\n";
} 

## Make output directory (and path) if it doesn't exist
my $orig_umask = umask;
umask 0000;
unless( -e $output or make_path($output , { verbose => 1, mode => 0775, } )) {
    umask $orig_umask;
    die "Unable to create $output\n";
}
umask $orig_umask;

my $FACETS = '';
my $PYTHON = '';
my $PERL = '';
my $VEP = '';
my $HG19_FASTA = '';
my $VCF2MAF = '';
my $B37_FASTA = '';
my $B37_MM10_HYBRID_FASTA = '';
my $MM10_FASTA = '';
my $MM9_FASTA = ''; 
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
    elsif($conf[0] =~ /facets/i){
        if(!-e "$conf[1]/facets"){
            die "CAN'T FIND facets IN $conf[1] $!";
        }
        $FACETS = $conf[1];
    }
    elsif($conf[0] =~ /^perl/i){
	if(!-e "$conf[1]/perl"){
	    die "CAN'T FIND perl IN $conf[1] $!";
	}
	$PERL = $conf[1];
    }
    elsif($conf[0] =~ /^r$/i){
        if(!-e "$conf[1]/R"){
            die "CAN'T FIND R IN $conf[1] $!";
        }
        my $path_tmp = $ENV{'PATH'};
        $ENV{'PATH'} = "$conf[1]:$path_tmp";
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
    elsif($conf[0] =~ /mm10_fasta/i){
        if(!-e "$conf[1]"){
            if($species =~ /mouse|mm10/i){
                die "CAN'T FIND $conf[1] $!";
            }
        }
        $MM10_FASTA = $conf[1];
    }
    elsif($conf[0] =~ /mm9_fasta/i){
        if(!-e "$conf[1]"){
            if($species =~ /mm9/i){
                die "CAN'T FIND $conf[1] $!";
            }
        }
        $MM9_FASTA = $conf[1];
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
     elsif($conf[0] =~ /hg19_fasta/i){
        if(!-e "$conf[1]"){
            die "CAN'T FIND $conf[1] $!";
        }
        $HG19_FASTA = $conf[1];
     }
}
close CONFIG;

# die if vcf2maf is empty, VEP, etc.
if(!-e $VCF2MAF ) {
    die "Cannot find vcf2maf. Either add or correct config file. $!";
}
if(!-e $VEP) {
    die "Cannot find vep. Either add or correct config file. $!";
} 

## FOR RIGHT NOW, ONLY HG19 IS BEING USED

my $REF_FASTA = $HG19_FASTA;
my $ncbi = "GRCh37";
if($species =~ /b37/i){
    $REF_FASTA = $B37_FASTA;
}elsif($species =~/mm10|mouse/i){
    $REF_FASTA = $MM10_FASTA;
    $ncbi = "GRCm38 --species mus_musculus";
}elsif($species =~/mm9/i){
    $REF_FASTA = $MM9_FASTA;
    $ncbi = "GRCm38 --species mus_musculus";
}
my $REF_DICT= $REF_FASTA;
$REF_DICT =~ s/\.[^.]+$//;
$REF_DICT = "$REF_DICT.dict";

##
## Softlink necessary files to output dir. Change the variables to the softlinked area.
##
my $base_hc = basename($hc_vcf);
if ( -l "$output/$base_hc" ) {
    unlink "$output/$base_hc";
}
symlink("$hc_vcf", "$output/$base_hc");

$hc_vcf = "$output/$base_hc";

for my $vcf (@mutect_vcfs){
    my $mutect_vcf = "$mutect_dir/$vcf";
    if ( -l "$output/" . basename($mutect_vcf) ) {
        unlink "$output/" . basename($mutect_vcf);
    }

    my $mutect_txt = $mutect_vcf;
    $mutect_txt =~ s/vcf$/txt/g;

    if ( -l "$output/" . basename($mutect_txt) ) {
        unlink "$output/" . basename($mutect_txt);
    }
    symlink($mutect_txt, "$output/" . basename($mutect_txt));
    symlink($mutect_vcf, "$output/" . basename($mutect_vcf));
}

my $progress = "$output/progress";
if (! -d "$progress"){
    mkdir("$progress", 0755) or die "Making progress directory did not work. $!";
}


# Beginning of commands

if($force_run || ! -e "$progress/" . basename("$hc_vcf") . "_vcf2maf0.done"){
    $force_run = 1;
    ##
    ##
    ## Get Indels from Haplotype Caller
    ##
    ## vcf2maf - Change HC output to maf format
    if($pairing){
        print "$PYTHON/python $Bin/maf/vcf2maf0.py -i $hc_vcf -c haplotypecaller -o $hc_vcf\_HC.maf0 -p $pairing\n\n";
        `$PYTHON/python $Bin/maf/vcf2maf0.py -i $hc_vcf -c haplotypecaller -o $hc_vcf\_HC.maf0 -p $pairing`;
    }
    else{
        print "$PYTHON/python $Bin/maf/vcf2maf0.py -i $hc_vcf -c haplotypecaller -o $hc_vcf\_HC.maf0\n\n";
        `$PYTHON/python $Bin/maf/vcf2maf0.py -i $hc_vcf -c haplotypecaller -o $hc_vcf\_HC.maf0`;
    }
    &checkResult($?, $progress,  basename("$hc_vcf") . "_vcf2maf0");
}

if($force_run || ! -e "$progress/" . basename("$hc_vcf") . "_pAqSomHC.done"){
    $force_run = 1;
    ## Quality Filtering for haplotype caller 
    print "$PYTHON/python $Bin/maf/pA_qSomHC.py < $hc_vcf\_HC.maf0 > $hc_vcf\_HC.maf1\n\n";
    `$PYTHON/python $Bin/maf/pA_qSomHC.py < $hc_vcf\_HC.maf0 > $hc_vcf\_HC.maf1`;
    &checkResult($?, $progress,  basename("$hc_vcf") . "_pAqSomHC");
}

if($force_run || ! -e "$progress/" . basename("$hc_vcf") . "_oldMaf2tcgaMaf.done"){
    $force_run = 1;
    ## Convert to TCGA MAF format
    print "$PYTHON/python $Bin/maf/oldMAF2tcgaMAF.py $species $hc_vcf\_HC.maf1 $hc_vcf\_HC.maf2\n\n";
    `$PYTHON/python $Bin/maf/oldMAF2tcgaMAF.py $species $hc_vcf\_HC.maf1 $hc_vcf\_HC.maf2`;
    &checkResult($?, $progress,  basename("$hc_vcf") . "_oldMaf2tcgaMaf");
}

if($force_run || ! -e "$progress/" . basename("$hc_vcf") . "_indelOnly.done"){
    $force_run = 1;
    ## Select for Indels Only
    print "$PYTHON/python $Bin/maf/indelOnly.py < $hc_vcf\_HC.maf2 > $hc_vcf\_qSomHC_InDels_TCGA_MAF.txt\n\n";
    `$PYTHON/python $Bin/maf/indelOnly.py < $hc_vcf\_HC.maf2 > $hc_vcf\_qSomHC_InDels_TCGA_MAF.txt`;
    &checkResult($?, $progress,  basename("$hc_vcf") . "_indelOnly");
}


##
##
## Get DMP re-filtered MAF from MuTect
##

foreach my $line (@pair_lines) {
    my($normal,$tumor) = split(/\s+/, $line);
    next if ($tumor eq 'na' || $normal eq 'na');
    my @vcfIndex = indices("$normal\_$tumor\_mutect", @mutect_vcfs);
    
    my $vcf = $mutect_vcfs[$vcfIndex[0]];
    my $linked_vcf = "$output/$vcf";
    my $linked_txt = $linked_vcf;
    $linked_txt =~ s/vcf$/txt/g;

    if($force_run || ! -e "$progress/" . basename("$vcf") . "_vcf2maf0.done"){
        $force_run = 1;
        ## Change vcf to maf, using the extra .txt file
        print "$PYTHON/python $Bin/maf/vcf2maf0.py -i $linked_vcf -aF $linked_txt -c mutect -o $linked_vcf\_maf0.txt -p $pairing -n $normal -t $tumor\n\n";
        `$PYTHON/python $Bin/maf/vcf2maf0.py -i $linked_vcf -aF $linked_txt -c mutect -o $linked_vcf\_maf0.txt -p $pairing -n $normal -t $tumor`;
        &checkResult($?, $progress,  basename("$vcf") . "_vcf2maf0"); 
    }

    if($force_run || ! -e "$progress/" . basename("$vcf") . "_DMP_rescue.done"){
        $force_run = 1;
        ## rescue 
        print "$PYTHON/python $Bin/rescue/DMP_rescue.py  < $linked_vcf\_maf0.txt > $linked_vcf\_maf1.txt 2> $linked_vcf\_maf_rescue.log\n\n";
        `$PYTHON/python $Bin/rescue/DMP_rescue.py  < $linked_vcf\_maf0.txt > $linked_vcf\_maf1.txt 2> $linked_vcf\_maf_rescue.log `;
        &checkResult($?, $progress,  basename("$vcf") . "_DMP_rescue");
    }

    if($force_run || ! -e "$progress/" . basename("$vcf") . "_oldMaf2tcgaMaf.done"){
        $force_run = 1;
        ## Change to tcga maf
        print "$PYTHON/python $Bin/maf/oldMAF2tcgaMAF.py $species $linked_vcf\_maf1.txt $linked_vcf\_maf2.txt\n\n";
        `$PYTHON/python $Bin/maf/oldMAF2tcgaMAF.py $species $linked_vcf\_maf1.txt $linked_vcf\_maf2.txt`;
        &checkResult($?, $progress,  basename("$vcf") . "_oldMaf2tcgaMaf");   
    

        ## Select for passing variants
        print "awk -F\"\\t\" '\$40==\"FILTER\"||\$40==\"PASS\"{print \$0}' $linked_vcf\_maf2.txt > $linked_vcf\_DMPFilter_TCGA_MAF.txt\n\n";
        `awk -F"\t" '\$40=="FILTER"||\$40=="PASS"{print \$0}' $linked_vcf\_maf2.txt > $linked_vcf\_DMPFilter_TCGA_MAF.txt`;
        &checkResult($?, $progress,  basename("$vcf") . "_oldMaf2tcgaMaf");
    } 
}

##
##
## Merge and Annotate?
##

if($force_run || ! -e "$progress/$pre\_merge_maf.done"){
    $force_run = 1;
    ## Get the first 39 fields of HC 
    print "cut -f-39 $hc_vcf\_qSomHC_InDels_TCGA_MAF.txt > $output/merge_maf0.txt\n\n";
    `cut -f-39 $hc_vcf\_qSomHC_InDels_TCGA_MAF.txt > $output/$pre\_merge_maf0.txt`;
    &checkResult($?, $progress, "$pre\_merge_maf");

    ## For each mutect result, get all results except for the header
    for my $vcf (@mutect_vcfs){
        my $linked_vcf = "$output/$vcf";
        print "tail -n +2 $linked_vcf\_DMPFilter_TCGA_MAF.txt >> $output/$pre\_merge_maf0.txt\n\n";
        `tail -n +2 $linked_vcf\_DMPFilter_TCGA_MAF.txt >> $output/$pre\_merge_maf0.txt`;
        &checkResult($?, $progress, "$pre\_merge_maf");
    }
} 

push @maf_header, "#version 2.4";
push @maf_header, "#SVN Revision: $svnRev";

if($species !~ /hg19|b37|mm10|mouse|human/i) { ###   |mm10|mouse/i) { uncomment later!
    sleep(2);
    print "cut -f-34 $output/$pre\_merge_maf0.txt > $output/$pre\_haplotect_TCGA_MAF.txt\n\n";
    `cut -f-34 $output/$pre\_merge_maf0.txt > $output/$pre\_haplotect_TCGA_MAF.txt`;
    print "Must be human or mouse(mm10 )species to continue.\n";
    &addHeader;
    if($delete_temp){
        &cleanUp;
    }
    exit 0;
}

print "\n#######\n#######\nStarting VEP. \n";
if($force_run || ! -e "$progress/$pre\_run_vep.done"){
    $force_run = 1;
    # these are names needed for the "retain-cols" option in VEP
    my $VEP_COLUMN_NAMES = "Center,Verification_Status,Validation_Status,Mutation_Status,Sequencing_Phase,Sequence_Source,Validation_Method,Score,BAM_file,Sequencer,Tumor_Sample_UUID,Matched_Norm_Sample_UUID,Caller";

    ## Create tmp and ref directory. Delete these later*
    if( ! -d "$output/tmp/" ){
        print "$output/tmp/ does not exist. Will create it now\n";
        mkdir("$output/tmp", 0755) or die "Making tmp didn't work $!";
    }
    if( ! -d "$output/ref/" ){
        print "$output/ref/ does not exist. Will create it now\n";
        mkdir("$output/ref", 0755) or die "Making tmp didn't work $!";
    }

    my $ref_base = basename($REF_FASTA);

    # softlink reference
    symlink($REF_FASTA, "$output/ref/$ref_base");
    symlink("$REF_FASTA.fai", "$output/ref/$ref_base.fai");

    ## vep-forks is at 4 because that is how many CPUs we ask for
    print "\n$PERL/perl $VCF2MAF/maf2maf.pl --tmp-dir $output/tmp --ref-fasta $output/ref/$ref_base --ncbi-build $ncbi --vep-forks 4 --vep-path $VEP --vep-data $VEP --retain-cols $VEP_COLUMN_NAMES --input-maf $output/$pre\_merge_maf0.txt --output-maf $output/$pre\_merge_maf0.VEP\n\n";

    `$PERL/perl $VCF2MAF/maf2maf.pl --tmp-dir $output/tmp --ref-fasta $output/ref/$ref_base --ncbi-build $ncbi --vep-forks 4 --vep-path $VEP --vep-data $VEP --retain-cols $VEP_COLUMN_NAMES --input-maf $output/$pre\_merge_maf0.txt --output-maf $output/$pre\_merge_maf0.VEP > $output/$pre\_merge_maf0.log 2>&1`;
    &checkResult($?, $progress, "$pre\_run_vep");
}
push @maf_header, "#VCF2MAF maf2maf.pl\tVEP: $VEP\tNCBI: $ncbi\tLOCATION: $VCF2MAF";


if($force_run || ! -e "$progress/$pre\_tcgaMaf_portalMaf.done"){
    $force_run = 1;
    #This removes any records that don't have a gene name at the front
    #`grep -v ^Unknown $output/$pre\_merge_maf0.VEP > $output/$pre\_merge_maf0.VEP`;
    #&checkResult($?, $progress,  basename("$vcf") . "_tcgaMaf_portalMaf");
    `cut -f-34 $output/$pre\_merge_maf0.VEP > $output/$pre\_haplotect_TCGA_MAF.txt`;
    &checkResult($?, $progress, "$pre\_tcgaMaf_portalMaf");

    ## creating MAF for cbio portal submission
    print "$PYTHON/python $Bin/maf/pA_reSortCols.py -i $output/$pre\_merge_maf0.VEP -f $Bin/maf/finalCols_PORTAL.txt -o $output/$pre\_haplotect_TCGA_PORTAL_MAF.txt\n\n";
    `$PYTHON/python $Bin/maf/pA_reSortCols.py -i $output/$pre\_merge_maf0.VEP -f $Bin/maf/finalCols_PORTAL.txt -o $output/$pre\_haplotect_TCGA_PORTAL_MAF.txt`;
    &checkResult($?, $progress, "$pre\_tcgaMaf_portalMaf");
}

if($patient && $bam_dir){
    if($force_run || ! -e "$progress/$pre\_maf_fillout.done"){
        $force_run = 1;
        print "Starting maf fillout\n";
        ## @bamList was populated above
        my $bam_inputs = join(" ", @bamList);

        print "$Bin/maf/fillout/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta $REF_FASTA $bam_inputs --output $output/$pre/_haplotect_TCGA_basecounts.txt --maf $output/$pre\_haplotect_TCGA_PORTAL_MAF.txt --filter_improper_pair 0\n\n";
        `$Bin/maf/fillout/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta $REF_FASTA $bam_inputs --output $output/$pre\_haplotect_TCGA_basecounts.txt --maf $output/$pre\_haplotect_TCGA_PORTAL_MAF.txt --filter_improper_pair 0`;
        &checkResult($?, $progress, "$pre\_maf_fillout");

        print "$PYTHON/python $Bin/maf/fillout/dmp2portalMAF -s $species -m $output/$pre\_haplotect_TCGA_PORTAL_MAF.txt -p $pairing -P $patient -c haplotect -b $output/$pre\_haplotect_TCGA_basecounts.txt -o $output/$pre\_haplotect_TCGA_PORTAL_MAF_fillout.txt\n\n";
        `$PYTHON/python $Bin/maf/fillout/dmp2portalMAF -s $species -m $output/$pre\_haplotect_TCGA_PORTAL_MAF.txt -p $pairing -P $patient -c haplotect -b $output/$pre\_haplotect_TCGA_basecounts.txt -o $output/$pre\_haplotect_TCGA_PORTAL_MAF_fillout.txt`;
        &checkResult($?, $progress, "$pre\_maf_fillout");
    }
}


###
##
## BEGINNING OF EXTRA COLUMNS ** This will be updated with a better way of 
## grabbing and adding this information, but for right now we will just add
## what Nick's script does
##
###

## run tri nulceotide script

if($force_run || ! -e "$progress/$pre\_bedtools_anno.done"){
    $force_run = 1;
    my $extraStuff = '';
    if ($species =~ /hg19|human|b37/i){
        $extraStuff =  " --target $Bin/targets/IMPACT410_$species/IMPACT410_$species\_targets_plus5bp.bed --targetname IMPACT_410";
    }
    print "$PERL/perl $Bin/maf/bedtools_annotations.pl --in_maf $output/$pre\_merge_maf0.VEP --species $species --output $output --config $config --fastq $extraStuff \n\n";
    `$PERL/perl $Bin/maf/bedtools_annotations.pl --in_maf $output/$pre\_merge_maf0.VEP --species $species --output $output --config $config --fastq $extraStuff`;
    &checkResult($?, $progress, "$pre\_bedtools_anno");
}

# exac annotate 
if($species =~ /hg19|b37|human/i){
    if($force_run || ! -e "$progress/$pre\_exac_anno.done"){
        $force_run = 1; 
        print "perl $Bin/maf/exact_annotate.pl --in_maf $output/$pre\_merge_maf0.VEP --species $species --output $output --config $config --data $Bin/data\n\n";
        `$PERL/perl $Bin/maf/exact_annotate.pl --in_maf $output/$pre\_merge_maf0.VEP --species $species --output $output --config $config --data $Bin/data`;
        &checkResult($?, $progress, "$pre\_exac_anno");
    }

##
## Merging of extra columns ** I haven't created a way to merge the columns while creating them, so I still
## have to use Nick's mkTaylorMAF.py, which I am going to rename to mergeExtraCols.py
##
    if($force_run || ! -e "$progress/$pre\_merge_cols.done"){
        $force_run = 1;
        print "$PYTHON/python $Bin/maf/mergeExtraCols.py --triNuc $output/TriNuc.txt --targeted $output/maf_targets.IMPACT_410 --exac $output/exact.vcf --maf $output/$pre\_merge_maf0.VEP > $output/$pre\_haplotect_VEP_MAF.txt\n\n";
        `$PYTHON/python $Bin/maf/mergeExtraCols.py --triNuc $output/TriNuc.txt --targeted $output/maf_targets.IMPACT_410 --exac $output/exact.vcf --maf $output/$pre\_merge_maf0.VEP > $output/$pre\_haplotect_VEP_MAF.txt`;
        &checkResult($?, $progress, "$pre\_merge_cols");
    }

} else {
    if($force_run || ! -e "$progress/$pre\_merge_cols.done"){
        $force_run = 1;
        `/bin/touch $output/blank`;
        print "$PYTHON/python $Bin/maf/mergeExtraCols.py --triNuc $output/TriNuc.txt --maf $output/$pre\_merge_maf0.VEP > $output/$pre\_haplotect_VEP_MAF.txt\n\n";
        `$PYTHON/python $Bin/maf/mergeExtraCols.py --triNuc $output/TriNuc.txt --maf $output/$pre\_merge_maf0.VEP > $output/$pre\_haplotect_VEP_MAF.txt`;
        &checkResult($?, $progress, "$pre\_merge_cols");
    }
}


sub checkResult{
    my ($status, $progress, $filebase) = @_;

    if($status == 0)
    {
        `/bin/touch $progress/$filebase.done`;
    } else{
        `/bin/rm -f $progress/$filebase.done`;
        print "\nThere was an error with script $filebase!\n";
        exit 1
    }

}

##
## Add headers to files that will be kept. 
##
## Print maf headers
open ( my $fh, '>', "$output/header.txt");
while(@maf_header){
    print $fh shift @maf_header;
    print $fh "\n";
}
close $fh;

##
## Add header
##
&addHeader;

##
## Delete Temp Files
##
if($delete_temp){
    &cleanUp;
}

sub addHeader{
    my @outFileArray = ();
    if($species !~ /hg19|b37|mm10|mouse|human/i) {
        @outFileArray = ( "$output/$pre\_haplotect_TCGA_MAF.txt");
    } elsif($patient && $bam_dir){
        @outFileArray = ( "$output/$pre\_haplotect_VEP_MAF.txt", "$output/$pre\_haplotect_TCGA_PORTAL_MAF_fillout.txt", "$output/$pre\_haplotect_TCGA_PORTAL_MAF.txt", "$output/$pre\_haplotect_TCGA_MAF.txt");
    } else {
        @outFileArray = ( "$output/$pre\_haplotect_VEP_MAF.txt", "$output/$pre\_haplotect_TCGA_PORTAL_MAF.txt", "$output/$pre\_haplotect_TCGA_MAF.txt" );
    }

    foreach my $f (@outFileArray){
        `grep -v ^# $f > $f\_temp`;
        `cat $output/header.txt $f\_temp > $f`;
        unlink("$f\_temp");
    }

}

sub cleanUp{
    my @files = glob( $output . '/*.maf[12]');
    for my $fname (@files) {
        print "Removing: $fname \n";
        unlink($fname);
    }
    @files = glob ("$output/exact.vcf $output/header.txt $output/maf_targets.* $output/TriNuc.txt $output/*mutect_calls* $output/blank $output/*Haplotype* $output/*maf2.txt* $output/ref $output/tmp $output/*basecounts.txt $output/$pre\_merge_maf0.VEP");


    for my $fname (@files) {
        print "Removing: $fname \n";
        unlink($fname);
    }
    `rm -r $output/xtra $output/bed $output/ref $output/tmp $progress`;
}

sub indices {
  my ($keyword, @paths) = @_;
  my @results = ();

  for ( my $i = 0; $i < @paths; ++$i )
  {
    if ($paths[$i] =~ /$keyword/)
    {
      push @results, $i;
    }
  }

  return @results;
}

