#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin);
use File::Path qw(make_path);
use File::Basename;

my ($svnRev, $pairing, $hc_vcf, $help, $mutect_dir, $species, $bam_dir, $patient, $pre, $output, $delete_temp, $config, $exac_vcf);

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
            'exac_vcf=s' => \$exac_vcf,
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
if($species !~ /hg19|b37|mm10|hybrid/i){
    print "Haplotect only works 100% with hg19/b37/mm10 for right now, otherwise we will print out the MAF before annotation and exit\n\n";
}

## Checking HC VCF
if( ! -e $hc_vcf || -z $hc_vcf){
    die "The haplotype vcf file either does not exist or is empty. ($hc_vcf)\n";
}

if($hc_vcf !~ /^\//){
    $hc_vcf = "$curDir/$hc_vcf";
}

## Checking exac vcf
if($exac_vcf){
    if( !-e $exac_vcf){
        die "$exac_vcf DOES NOT EXIST";
    }
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
my $FIXMULTIINDEL = '';
my $VCFTOOLS = '';
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
    } elsif($conf[0] =~ /facets/i){
        if(!-e "$conf[1]/facets"){
            die "CAN'T FIND facets IN $conf[1] $!";
        }
        $FACETS = $conf[1];
    }elsif($conf[0] =~/fixmultiindel/i){
        if(!-e "$conf[1]/fixMultiInDel.sh"){
            die "CAN'T FIND fixMultiInDel in $conf[1] $!";
        }
        $FIXMULTIINDEL = $conf[1];
    } elsif($conf[0] =~ /^perl/i){
        if(!-e "$conf[1]/perl"){
            die "CAN'T FIND perl IN $conf[1] $!";
        }
	    $PERL = $conf[1];
    } elsif($conf[0] =~ /^r$/i){
        if(!-e "$conf[1]/R"){
            die "CAN'T FIND R IN $conf[1] $!";
        }
        my $path_tmp = $ENV{'PATH'};
        $ENV{'PATH'} = "$conf[1]:$path_tmp";
    } elsif($conf[0] =~ /b37_fasta/i){
        if(!-e "$conf[1]"){
            if($species =~ /human|^b37$/i){
                die "CAN'T FIND $conf[1] $!";
            }
        }
        $B37_FASTA = $conf[1];
    } elsif($conf[0] =~ /b37_mm10_hybrid_fasta/i){
        if(!-e "$conf[1]"){
            if($species =~ /hybrid|b37_mm10/i){
                die "CAN'T FIND $conf[1] $!";
            }
        }
        $B37_MM10_HYBRID_FASTA = $conf[1];
    } elsif($conf[0] =~ /mm10_fasta/i){
        if(!-e "$conf[1]"){
            if($species =~ /mouse|mm10/i){
                die "CAN'T FIND $conf[1] $!";
            }
        }
        $MM10_FASTA = $conf[1];
    } elsif($conf[0] =~ /mm9_fasta/i){
        if(!-e "$conf[1]"){
            if($species =~ /mm9/i){
                die "CAN'T FIND $conf[1] $!";
            }
        }
        $MM9_FASTA = $conf[1];
    } elsif($conf[0] =~ /^vep/i){
        if(!-e "$conf[1]/variant_effect_predictor.pl"){
            die "CAN'T FIND VEP IN $conf[1] $!";
        }
        $VEP = $conf[1];
    } elsif($conf[0] =~/vcf2maf/i){
        if(!-e "$conf[1]/vcf2maf.pl"){
            die "CAN'T FIND vcf2maf.pl in $conf[1] $!";
        }
        $VCF2MAF = $conf[1];
    } elsif($conf[0] =~/vcftools/i){
        if(!-e "$conf[1]/vcftools"){
            die "CAN'T FIND vcftools in $conf[1] $!";
        }
        $VCFTOOLS = $conf[1];
    } elsif($conf[0] =~ /bcftools/i){
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
    } elsif($conf[0] =~ /samtools/i){
        if(!-e "$conf[1]/samtools"){
            die "CAN'T FIND samtools IN $conf[1] $!";
        }
        my $path_tmp = $ENV{'PATH'};
        $ENV{'PATH'} = "$conf[1]:$path_tmp";
     } elsif($conf[0] =~ /tabix/i){
        if(!-e "$conf[1]/tabix"){
            die "CAN'T FIND tabix IN $conf[1] $!";
        }
        my $path_tmp = $ENV{'PATH'};
        $ENV{'PATH'} = "$conf[1]:$path_tmp";
     } elsif($conf[0] =~ /hg19_fasta/i){
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
}elsif($species =~ /hybrid/){
    $species = "b37";
    $REF_FASTA = $B37_MM10_HYBRID_FASTA;
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
# Do the paired part of the analysis from the generate maf for haplotect. 
# Then follow up and do the mutect stuff too.

# To use vcf2maf we first have to split the vcf to pairs. This will be output in $output/tmp
if( ! -d "$output/tmp/" ){
    print "$output/tmp/ does not exist. Will create it now\n";
    mkdir("$output/tmp", 0755) or die "Making tmp dir didn't work $!";
}
    
my @indv_mafs;
my @mutect_mafs;

## Haplotype Caller pre processing. Split and select for indels!!
print "$FIXMULTIINDEL/fixMultiInDel.sh  $hc_vcf $hc_vcf\_split.vcf";
`$FIXMULTIINDEL/fixMultiInDel.sh  $hc_vcf $hc_vcf\_split.vcf`;

$hc_vcf = "$hc_vcf\_split.vcf";

foreach my $line (@pair_lines) {
    my($normal,$tumor) = split(/\s+/, $line);
    next if ($tumor eq 'na' || $normal eq 'na');

    ###
    ### HAPLOTYPE CALLER PART 1 (pair wise running)
    ###

    ## split vcf into paired vcfs     
    my $split_vcf_prefix = "$output/tmp/" . basename($hc_vcf) . "_$tumor\_$normal";  
    if($force_run || !-e "$progress/" . basename("$hc_vcf") . "_$tumor\_$normal\_vcf.done"){
        $force_run = 1;
        &splitVcf($hc_vcf, $split_vcf_prefix, $tumor, $normal, "$progress/",  basename("$hc_vcf") . "_$tumor\_$normal\_vcf", $split_vcf_prefix . ".recode.vcf");    
    }

    ## Filtering for actual variants (since we are splitting a vcf and the vcftools puts all records in vcf output)
    if($force_run || !-e "$progress/" . basename("$hc_vcf") . "_$tumor\_$normal\_vcf_filterNonVars.done"){
        $force_run = 1;
        print "$PYTHON/python $Bin/maf/filter_nonVars_vcf.py -v $split_vcf_prefix\.recode.vcf -t $tumor -n $normal --somatic -o $split_vcf_prefix\_vars.vcf \n";
        `$PYTHON/python $Bin/maf/filter_nonVars_vcf.py -v $split_vcf_prefix\.recode.vcf -t $tumor -n $normal --somatic -o $split_vcf_prefix\_vars.vcf`;
        &checkResult($?, $progress, basename("$hc_vcf") . "_$tumor\_$normal\_vcf_filterNonVars")
    }
    ## run vcf2maf
    if($force_run || !-e "$progress/" . basename("$hc_vcf") . "_$tumor\_$normal\_vcf2maf.done"){
        $force_run = 1;
        &run_vcf_2_maf("$split_vcf_prefix\_vars.vcf", $split_vcf_prefix . ".maf", $tumor, $normal,  basename("$hc_vcf") . "_$tumor\_$normal\_vcf2maf");
    }

    # TODO: Instead of just assuming the correct column is $10, find it first, then use it
    
    # should this be in the regular folder, or in the tmp folder
    my $indel_maf = $split_vcf_prefix . "_indel_MAF.txt";
    if($force_run || !-e "$progress/" . basename("$split_vcf_prefix") . "_combine_MAF.done"){
        $force_run = 1;
         print "grep -e ^# -e ^Hugo_Symbol $split_vcf_prefix\.maf > $indel_maf\n";
         `grep -e ^# -e ^Hugo_Symbol $split_vcf_prefix\.maf > $indel_maf`;
         `awk 'BEGIN{FS="\t"} \$10=="INS" || \$10=="DEL" ' $split_vcf_prefix\.maf >> $indel_maf `;
         
         &checkResult($?, $progress,  basename("$split_vcf_prefix") . "_combine_MAF");
           
    }
    push(@indv_mafs, $indel_maf); 
        
    ###
    ### MUTECT PART
    ###
    
    my @vcfIndex = indices("$normal\_$tumor\_mutect", @mutect_vcfs);
    
    my $vcf = $mutect_vcfs[$vcfIndex[0]];
    my $linked_vcf = "$output/$vcf";
    my $linked_txt = $linked_vcf;
    $linked_txt =~ s/vcf$/txt/g;
    
    # Do resuce of failed mutect stuff

    if($force_run || ! -e "$progress/$vcf\_DMP_rescue.done"){
        $force_run = 1;
        ## rescue 
        # TODO: Does this script add "caller"? If not then it should!
        print "$PYTHON/python $Bin/rescue/DMP_rescue.py  -v $linked_vcf --txt $linked_txt -t $tumor -n $normal -o $output/tmp/$vcf\_rescued.vcf 2> $output/tmp/$vcf\_maf_rescue.log\n\n";
        `$PYTHON/python $Bin/rescue/DMP_rescue.py  -v $linked_vcf --txt $linked_txt -t $tumor -n $normal -o $output/tmp/$vcf\_rescued.vcf 2> $output/tmp/$vcf\_maf_rescue.log `;
        &checkResult($?, $progress,  "$vcf\_DMP_rescue");
    }
    
    # vcf2maf
    ## run vcf2maf
    if($force_run || !-e "$progress/$vcf\_$tumor\_$normal\_mutect_vcf2maf.done"){
        $force_run = 1;
        &run_vcf_2_maf("$output/tmp/$vcf\_rescued.vcf",  "$output/tmp/$vcf\_mutect.maf", $tumor, $normal,  "$vcf\_$tumor\_$normal\_mutect_vcf2maf");
    }
    
    # TODO: The same as before, make sure that $40 is really "FILTER", always.
    
    # now only save passing variants (pass or rescue)
    my $passing_maf = "$output/tmp/$vcf\_mutect_passed.maf";
    if($force_run || !-e "$progress/$vcf\_$tumor\_$normal\_mutect_rm_failed.done"){
    	`grep -e ^# -e ^Hugo_Symbol $output/tmp/$vcf\_mutect.maf > $passing_maf`;
        my $filter_index = `grep ^Hugo_Symbol $output/tmp/$vcf\_mutect.maf | tr "\t" "\n" | grep -n "^FILTER\$" | cut -d ":" -f 1`;
        chomp $filter_index;
        if($filter_index eq ""){
            die "Cannot find FILTER column in maf $output/tmp/$vcf\_mutect.maf $!";
        }
        ` awk 'BEGIN{FS="\t"} \$$filter_index == "PASS" || \$$filter_index == "RESCUE" ' $output/tmp/$vcf\_mutect.maf >> $passing_maf\_temp `;
        ` awk 'BEGIN{FS="\t"} {if(\$$filter_index == "FILTER") { print \$0"\tCaller" } else if (\$$filter_index == "RESCUE") {print \$0"\tmutect.Rescue"} else if(\$$filter_index == "PASS") {print \$0"\tmutect"} else {print \$0}} ' $passing_maf\_temp > $passing_maf`;
    	&checkResult($?, $progress,  "$vcf\_$tumor\_$normal\_mutect_rm_failed");
    }
    push(@mutect_mafs, $passing_maf);
}

###
### Pairs are done! Now Merge both Haplotype caller together (then filter) and combine that and mutect
###

#COMBINE SNP FILES
my $mut_maf = "$output/$pre\_mutect_SNPs.maf";
if($force_run || !-e "$progress/$pre\_combine_SNPS.done"){
    $force_run = 1;
    my $mafs = join(" ", @mutect_mafs);
    print "grep -e ^# -e ^Hugo_Symbol $indv_mafs[0] > $mut_maf\n";
    print "cat $mafs | grep -v -e ^# -e ^Hugo_Symbol >> $mut_maf\n";
    `grep -e ^# -e ^Hugo_Symbol $indv_mafs[0] > $mut_maf`;
    `cat $mafs | grep -v -e ^# -e ^Hugo_Symbol >> $mut_maf`;
    
    &checkResult($?, $progress,  "$pre\_combine_SNPS");
}

#COMBINE INDEL FILES
my $raw_maf = "$output/$pre\_haplotect_RAW_INDEL.maf";
if($force_run || !-e "$progress/$pre\_combine_INDELS.done"){
    $force_run = 1;
    # Take array of all tumor/normal pairs, and add them together.
    my $mafs = join(" ", @indv_mafs);
    
    print 'grep -e ^# -e ^Hugo_Symbol $indv_mafs[0] > $raw_maf\n';
    print "cat $mafs | grep -v -e ^# -e ^Hugo_Symbol >> $raw_maf\n";
    `grep -e ^# -e ^Hugo_Symbol $indv_mafs[0] > $raw_maf`;
    `cat $mafs | grep -v -e ^# -e ^Hugo_Symbol >> $raw_maf`;

     &checkResult($?, $progress,  "$pre\_combine_INDELS");
}

# then filter INDELS based on somatic contraints
if($force_run || ! -e "$progress/$pre\_haplotect_pA_qFiltersHC.done"){
    $force_run=1;
        
    print "$PYTHON/python $Bin/maf/pA_qFiltersHC.py -m $raw_maf -s -c HaplotypeCaller -o $output/$pre\_haplotect_BIC_INDELS.maf";
    `$PYTHON/python $Bin/maf/pA_qFiltersHC.py -m $raw_maf -s -c HaplotypeCaller -o $output/$pre\_haplotect_BIC_INDELS.maf`;
    &checkResult($?, $progress, "$pre\_haplotect_pA_qFiltersHC");
}


##
## Merge and Annotate?
##

if($force_run || ! -e "$progress/$pre\_merge_maf.done"){
    $force_run = 1;
    ## Get the first 39 fields of HC 
    print "cp $output/$pre\_haplotect_BIC_INDELS.maf $output/merge_maf0.txt\n\n";
    `cp $output/$pre\_haplotect_BIC_INDELS.maf $output/$pre\_merge_maf0.txt`;
    &checkResult($?, $progress, "$pre\_merge_maf");

    # Now get the mutect
    print "grep -v -e ^# -e ^Hugo_Symbol $mut_maf >> $output/$pre\_merge_maf0.txt\n\n";
    `grep -v -e ^# -e ^Hugo_Symbol $mut_maf >> $output/$pre\_merge_maf0.txt`;
    &checkResult($?, $progress, "$pre\_merge_maf");

} 

push @maf_header, "#version 2.4";
push @maf_header, "#SVN Revision: $svnRev";

if($species !~ /hg19|b37|mm10|mouse|human|hybrid/i) { ###   |mm10|mouse/i) { uncomment later!
    sleep(2);
    print "$PYTHON/python $Bin/maf/pA_reSortCols.py -i $output/$pre\_merge_maf0.txt -f $Bin/maf/TCGA_MAF_fields.txt -o $output/$pre\_haplotect_TCGA_MAF.txt\n\n";
    `$PYTHON/python $Bin/maf/pA_reSortCols.py -i $output/$pre\_merge_maf0.txt -f $Bin/maf/TCGA_MAF_fields.txt -o $output/$pre\_haplotect_TCGA_MAF.txt`;
    print "Must be human or mouse(mm10 )species to continue.\n";
    &addHeader;
    if($delete_temp){
        &cleanUp;
    }
    exit 0;
}


if($force_run || ! -e "$progress/$pre\_tcgaMaf_portalMaf.done"){
    $force_run = 1;
    #This removes any records that don't have a gene name at the front
    #`grep -v ^Unknown $output/$pre\_merge_maf0.txt > $output/$pre\_merge_maf0.txt`;
    #&checkResult($?, $progress,  basename("$vcf") . "_tcgaMaf_portalMaf");
    print "$PYTHON/python $Bin/maf/pA_reSortCols.py -i $output/$pre\_merge_maf0.txt -f $Bin/maf/TCGA_MAF_fields.txt -o $output/$pre\_haplotect_TCGA_MAF.txt\n\n";
    `$PYTHON/python $Bin/maf/pA_reSortCols.py -i $output/$pre\_merge_maf0.txt -f $Bin/maf/TCGA_MAF_fields.txt -o $output/$pre\_haplotect_TCGA_MAF.txt`;
    &checkResult($?, $progress, "$pre\_tcgaMaf_portalMaf");

    ## creating MAF for cbio portal submission
    print "$PYTHON/python $Bin/maf/pA_reSortCols.py -i $output/$pre\_merge_maf0.txt -f $Bin/maf/finalCols_PORTAL.txt -o $output/$pre\_haplotect_TCGA_PORTAL_MAF.txt\n\n";
    `$PYTHON/python $Bin/maf/pA_reSortCols.py -i $output/$pre\_merge_maf0.txt -f $Bin/maf/finalCols_PORTAL.txt -o $output/$pre\_haplotect_TCGA_PORTAL_MAF.txt`;
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
    if ($species =~ /hg19|human|b37|hybrid/i){
        $extraStuff =  " --target $Bin/targets/IMPACT410_$species/IMPACT410_$species\_targets_plus5bp.bed --targetname IMPACT_410";
    }
    print "$PERL/perl $Bin/maf/bedtools_annotations.pl --in_maf $output/$pre\_merge_maf0.txt --species $species --output $output --config $config --fastq $extraStuff \n\n";
    `$PERL/perl $Bin/maf/bedtools_annotations.pl --in_maf $output/$pre\_merge_maf0.txt --species $species --output $output --config $config --fastq $extraStuff`;
    &checkResult($?, $progress, "$pre\_bedtools_anno");
}

# exac annotate 
if($species =~ /hg19|b37|human|hybrid/i){
    if($force_run || ! -e "$progress/$pre\_exac_anno.done"){
        $force_run = 1; 
        print "perl $Bin/maf/exac_annotate.pl --in_maf $output/$pre\_merge_maf0.txt --species $species --output $output --config $config --data $Bin/data\n\n";
        `$PERL/perl $Bin/maf/exac_annotate.pl --in_maf $output/$pre\_merge_maf0.txt --species $species --output $output --config $config --data $Bin/data`;
        &checkResult($?, $progress, "$pre\_exac_anno");
    }

##
## Merging of extra columns ** I haven't created a way to merge the columns while creating them, so I still
## have to use Nick's mkTaylorMAF.py, which I am going to rename to mergeExtraCols.py
##
    if($force_run || ! -e "$progress/$pre\_merge_cols.done"){
        $force_run = 1;
        print "$PYTHON/python $Bin/maf/mergeExtraCols.py --triNuc $output/TriNuc.txt --targeted $output/maf_targets.IMPACT_410 --exac $output/exac.vcf --maf $output/$pre\_merge_maf0.txt > $output/$pre\_haplotect_VEP_MAF.txt\n\n";
        `$PYTHON/python $Bin/maf/mergeExtraCols.py --triNuc $output/TriNuc.txt --targeted $output/maf_targets.IMPACT_410 --exac $output/exac.vcf --maf $output/$pre\_merge_maf0.txt > $output/$pre\_haplotect_VEP_MAF.txt`;
        &checkResult($?, $progress, "$pre\_merge_cols");
    }

} else {
    if($force_run || ! -e "$progress/$pre\_merge_cols.done"){
        $force_run = 1;
        `/bin/touch $output/blank`;
        print "$PYTHON/python $Bin/maf/mergeExtraCols.py --triNuc $output/TriNuc.txt --maf $output/$pre\_merge_maf0.txt > $output/$pre\_haplotect_VEP_MAF.txt\n\n";
        `$PYTHON/python $Bin/maf/mergeExtraCols.py --triNuc $output/TriNuc.txt --maf $output/$pre\_merge_maf0.txt > $output/$pre\_haplotect_VEP_MAF.txt`;
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
    @files = glob ("$output/exac.vcf $output/header.txt $output/maf_targets.* $output/TriNuc.txt $output/*mutect_calls* $output/blank $output/*Haplotype* $output/*maf2.txt* $output/ref $output/tmp $output/*basecounts.txt $output/$pre\_merge_maf0.txt");


    for my $fname (@files) {
        print "Removing: $fname \n";
        unlink($fname);
    }
    `rm -r $output/xtra $output/bed $output/ref $output/tmp $progress`;
}

sub splitVcf{
    my ($inVCF, $outPrefix, $tumor, $normal, $progress, $doneFile, $checkFile) = @_;
    #Now I have to use vcftools to separate the paired stuff
    print "$VCFTOOLS/vcftools --indv $tumor --indv $normal --vcf $inVCF --out $outPrefix --recode\n";
    `$VCFTOOLS/vcftools --indv $tumor --indv $normal --vcf $inVCF --out $outPrefix --recode`;
    &checkResult($?, $progress, $doneFile, $checkFile);
}

sub run_vcf_2_maf{
    my ($in_vcf, $out_maf, $tumor, $normal, $doneFile) = @_;
    
    if( ! -d "$output/ref/" ){
        print "$output/ref/ does not exist. Will create it now\n";
        mkdir("$output/ref", 0755) or die "Making ref didn't work $!";
    }
    my $ref_base = basename($REF_FASTA);
    
    # softlink reference
    symlink($REF_FASTA, "$output/ref/$ref_base");
    symlink("$REF_FASTA.fai", "$output/ref/$ref_base.fai");
    my $addOptions = '';
    if($exac_vcf){
        $addOptions = "--filter-vcf $exac_vcf";
    }
    
    print "\n#######\n#######\nStarting VEP. \n";
    print "$PERL/perl $VCF2MAF/vcf2maf.pl --input-vcf $in_vcf --species $species --output-maf $out_maf --ref-fasta $output/ref/$ref_base --tmp-dir $output/tmp/ --ncbi $ncbi --vep-forks 4 --vep-path $VEP --vep-data $VEP $addOptions --tumor-id $tumor --normal-id $normal \n";
    `$PERL/perl $VCF2MAF/vcf2maf.pl --input-vcf $in_vcf --species $species --output-maf $out_maf --ref-fasta $output/ref/$ref_base --tmp-dir $output/tmp/ --ncbi $ncbi --vep-forks 4 --vep-path $VEP --vep-data $VEP $addOptions --tumor-id $tumor --normal-id $normal`;

    &checkResult($?, $progress, $doneFile, $out_maf);
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



