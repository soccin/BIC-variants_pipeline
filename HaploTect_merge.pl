#!/opt/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin);
use File::Path qw(make_path);
use File::Basename;


my $ONCOTATORDIR="/opt/common/oncotator/oncotator_2014.01.06";

my ($pairing, $hc_vcf, $help, $mutect_dir, $pre, $output, $delete_temp);

# Default output = results/variants/FinalReport
# get current directory
$output = "results/variants/FinalReport";
my $curDir = `pwd`;
chomp $curDir;
my $commandline = join " ", "\n", $0, @ARGV;

# Print command line.
print "$commandline\n\n";

GetOptions ('pair=s' => \$pairing,
            'hc_vcf=s' => \$hc_vcf,
            'mutect_dir=s' => \$mutect_dir,
            'pre=s' => \$pre,
            'merge_dir=s' => \$output,
            'delete_temp' => \$delete_temp,
            'help|h' => \$help ) or exit(1);

# if the basic necessities are missing, exit
if(!$hc_vcf || !$mutect_dir || !$pre || $help){
    print <<HELP;
    
    USAGE: HaploTect_merge.pl -pair Proj_XXXX_sample_pairing.txt -hc_vcf Proj_XXXX_HaplotypeCaller.vcf -mutect_dir /ifs/work/mutect_results -pre Proj_XXXX -output /ifs/work/output 
    
    -pair\tFile listing tumor/normal pairing of samples for mutect/maf conversion; if not specified, considered unpaired
    -hc_vcf\tHaplotype caller's vcf output file
    -mutect_dir\tDirectory were MuTect output is found
    -pre\tProject ID Ex: Proj_4500
    -output\toutput directory where the merge intermediate and final files will be (default : results/variants/FinalReport)
    -delete_temp\tDelete the intermediate files.
    
HELP
    exit;
}

## Checking HC VCF
if( ! -e $hc_vcf || -z $hc_vcf){
    die "The haplotype vcf file either does not exist or is empty. ($hc_vcf)\n";
}
## Checking Mutect Directory
if( ! -d $mutect_dir ){
    die "The mutect directory either does not exist or is not a directory. ($mutect_dir)\n";
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

## Make output directory (and path) if it doesn't exist
my $orig_umask = umask;
umask 0000;
unless( -e $output or make_path($output , { verbose => 1, mode => 0775, } )) {
    umask $orig_umask;
    die "Unable to create $output\n";
}
umask $orig_umask;

##
## Softlink necessary files to output dir. Change the variables to the softlinked area.
##
my $base_hc = basename($hc_vcf);
if ( -l "$output/$base_hc" ) {
    unlink "$output/$base_hc";
}
symlink("$curDir/$hc_vcf", "$output/$base_hc");

$hc_vcf = "$output/$base_hc";

for my $vcf (@mutect_vcfs){
    my $mutect_vcf = "$curDir/$mutect_dir/$vcf";
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


##
##
## Get Indels from Haplotype Caller
##
## vcf2maf - Change HC output to maf format
if($pairing){
    print "$Bin/maf/vcf2maf0.py -i $hc_vcf -c haplotypecaller -o $hc_vcf\_HC.maf0 -p $pairing\n\n";
    `$Bin/maf/vcf2maf0.py -i $hc_vcf -c haplotypecaller -o $hc_vcf\_HC.maf0 -p $pairing`;
}
else{
    print "$Bin/maf/vcf2maf0.py -i $hc_vcf -c haplotypecaller -o $hc_vcf\_HC.maf0\n\n";
    `$Bin/maf/vcf2maf0.py -i $hc_vcf -c haplotypecaller -o $hc_vcf\_HC.maf0`;
}

## Quality Filtering for haplotype caller 
print "$Bin/maf/pA_qSomHC.py < $hc_vcf\_HC.maf0 > $hc_vcf\_HC.maf1\n\n";
`$Bin/maf/pA_qSomHC.py < $hc_vcf\_HC.maf0 > $hc_vcf\_HC.maf1`;

## Convert to TCGA MAF format
print "$Bin/maf/oldMAF2tcgaMAF.py hg19 $hc_vcf\_HC.maf1 $hc_vcf\_HC.maf2\n\n";
`$Bin/maf/oldMAF2tcgaMAF.py hg19 $hc_vcf\_HC.maf1 $hc_vcf\_HC.maf2`;

## Select for Indels Only
print "$Bin/maf/indelOnly.py < $hc_vcf\_HC.maf2 > $hc_vcf\_qSomHC_InDels_TCGA_MAF.txt\n\n";
`$Bin/maf/indelOnly.py < $hc_vcf\_HC.maf2 > $hc_vcf\_qSomHC_InDels_TCGA_MAF.txt`;



##
##
## Get DMP re-filtered MAF from MuTect
##

for my $vcf (@mutect_vcfs){
    my $linked_vcf = "$output/$vcf";
    my $linked_txt = $linked_vcf;
    $linked_txt =~ s/vcf$/txt/g;

    $_=$vcf;
    /_(s_.*?)_(s_.*?)_mutect/;
    my $normal="$1";
    my $tumor="$2";

    ## Change vcf to maf, using the extra .txt file
    print "$Bin/maf/vcf2maf0.py -i $linked_vcf -aF $linked_txt -c mutect -o $linked_vcf\_maf0.txt -p $pairing -n $normal -t $tumor\n\n";
    `$Bin/maf/vcf2maf0.py -i $linked_vcf -aF $linked_txt -c mutect -o $linked_vcf\_maf0.txt -p $pairing -n $normal -t $tumor`;
    
    ## rescue 
    print "$Bin/rescue/DMP_rescue.py  < $linked_vcf\_maf0.txt > $linked_vcf\_maf1.txt 2> $linked_vcf\_maf_rescue.log\n\n";
    `$Bin/rescue/DMP_rescue.py  < $linked_vcf\_maf0.txt > $linked_vcf\_maf1.txt 2> $linked_vcf\_maf_rescue.log `;
    
    ## Change to tcga maf
    print "$Bin/maf/oldMAF2tcgaMAF.py hg19 $linked_vcf\_maf1.txt $linked_vcf\_maf2.txt\n\n";
    `$Bin/maf/oldMAF2tcgaMAF.py hg19 $linked_vcf\_maf1.txt $linked_vcf\_maf2.txt`;
    
    ## Select for passing variants
    print "awk -F\"\\t\" '\$40==\"FILTER\"||\$40==\"PASS\"{print \$0}' $linked_vcf\_maf2.txt > $linked_vcf\_DMPFilter_TCGA_MAF.txt\n\n";
    `awk -F"\t" '\$40=="FILTER"||\$40=="PASS"{print \$0}' $linked_vcf\_maf2.txt > $linked_vcf\_DMPFilter_TCGA_MAF.txt`;
    
}


##
##
## Merge and Annotate?
##
## Get the first 39 fields of HC 
print "cut -f-39 $hc_vcf\_qSomHC_InDels_TCGA_MAF.txt > $output/merge_maf0.txt\n\n";
`cut -f-39 $hc_vcf\_qSomHC_InDels_TCGA_MAF.txt > $output/merge_maf0.txt`;

## For each mutect result, get all results except for the header
for my $vcf (@mutect_vcfs){
    my $linked_vcf = "$output/$vcf";
    print "tail -n +2 $linked_vcf\_DMPFilter_TCGA_MAF.txt >> $output/merge_maf0.txt\n\n";
    `tail -n +2 $linked_vcf\_DMPFilter_TCGA_MAF.txt >> $output/merge_maf0.txt`;
}

## link oncotator and oncotate!
symlink("$ONCOTATORDIR/db.properties", "./db.properties");
print "$ONCOTATORDIR/oncotateMaf.sh $output/merge_maf0.txt $output/merge_maf1.txt\n\n";
`$ONCOTATORDIR/oncotateMaf.sh $output/merge_maf0.txt $output/merge_maf1.txt`;

## fix Hugo?
print "$Bin/maf/pA_fixHugo.py <$output/merge_maf1.txt >$output/merge_maf2.txt\n\n";
`$Bin/maf/pA_fixHugo.py <$output/merge_maf1.txt >$output/merge_maf2.txt`;

## Update hugo symbol
print "$Bin/update_gene_names_and_ids.pl $output/merge_maf2.txt\n\n";
`/bin/rm $Bin/lib/hugo_data.tsv`;
`$Bin/update_gene_names_and_ids.pl $output/merge_maf2.txt > $output/merge_maf2_hugo.log 2>&1`;

## Grep for functional events YO!
print "cat $output/merge_maf2.txt_hugo_modified | $Bin/maf/pA_Functional_Oncotator2.py > $output/$pre\_HC_qSomHC_mutect_rescue_functional_MAF.txt\n\n";
`cat $output/merge_maf2.txt_hugo_modified | $Bin/maf/pA_Functional_Oncotator2.py > $output/$pre\_HC_qSomHC_mutect_rescue_functional_MAF.txt`;

## creating MAF for cbio portal submission
print "cat $output/$pre\_HC_qSomHC_mutect_rescue_functional_MAF.txt | $Bin/maf/pA_reSortCols.py $Bin/maf/finalCols_PORTAL.txt > $output/$pre\_HC_qSomHC_mutect_rescue_functional_MAF_PORTAL.txt\n\n";
`cat $output/$pre\_HC_qSomHC_mutect_rescue_functional_MAF.txt | $Bin/maf/pA_reSortCols.py $Bin/maf/finalCols_PORTAL.txt > $output/$pre\_HC_qSomHC_mutect_rescue_functional_MAF_PORTAL.txt`;


##
##
## Delete Temp Files
##
if($delete_temp){
    my @files = glob( $output . '/*maf[012].txt' );
    for my $fname (@files) {
        print "Removing: $fname \n";
        unlink($fname);
    }
}





