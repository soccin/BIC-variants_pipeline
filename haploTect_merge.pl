#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin);
use File::Path qw(make_path);
use File::Basename;

my ($pairing, $hc_vcf, $help, $mutect_dir, $bam_dir, $patient, $pre, $output, $delete_temp, $config);

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
            'align_dir=s' => \$bam_dir,
            'patient=s' => \$patient,
            'pre=s' => \$pre,
            'output=s' => \$output,
            'delete_temp' => \$delete_temp,
	    'config=s' => \$config,
            'help|h' => \$help ) or exit(1);

# if the basic necessities are missing, exit
if(!$hc_vcf || !$mutect_dir || !$pre || !$config || $help){
    print <<HELP;
    
    USAGE: HaploTect_merge.pl -pair Proj_XXXX_sample_pairing.txt -hc_vcf Proj_XXXX_HaplotypeCaller.vcf -mutect_dir /ifs/work/mutect_results -pre Proj_XXXX -output /ifs/work/output 
    
    -pair\tFile listing tumor/normal pairing of samples for mutect/maf conversion; if not specified, considered unpaired
    -hc_vcf\tHaplotype caller's vcf output file
    -mutect_dir\tDirectory were MuTect output is found
    -config\tConfiguration file
    -pre\tProject ID Ex: Proj_4500
    -patient\tPatient file (for optional fillout)
    -align_dir\tDirectory where the aligned bams are found (for optional fillout)
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

my $ONCOTATOR = '';
my $PYTHON = '';
my $PERL = '';
my $VEP = '';
my $HG19_FASTA = '';

open(CONFIG, "$config") or warn "CAN'T OPEN CONFIG FILE $config SO USING DEFAULT SETTINGS";
while(<CONFIG>){
    chomp;

    my @conf = split(/\s+/, $_);
    if($conf[0] =~ /oncotator/i){
	if(!-e "$conf[1]/oncotateMaf.sh"){
	    die "CAN'T FIND oncotateMaf.sh IN $conf[1] $!";
	}
	$ONCOTATOR = $conf[1];
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
    elsif($conf[0] =~ /hg19_fasta/i){
        if(!-e "$conf[1]"){
            die "CAN'T FIND $conf[1] $!";
        }
        $HG19_FASTA = $conf[1];
    }
}
close CONFIG;

## FOR RIGHT NOW, ONLY HG19 IS BEING USED

my $REF_FASTA = $HG19_FASTA;


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
    print "$PYTHON/python $Bin/maf/vcf2maf0.py -i $hc_vcf -c haplotypecaller -o $hc_vcf\_HC.maf0 -p $pairing\n\n";
    `$PYTHON/python $Bin/maf/vcf2maf0.py -i $hc_vcf -c haplotypecaller -o $hc_vcf\_HC.maf0 -p $pairing`;
}
else{
    print "$PYTHON/python $Bin/maf/vcf2maf0.py -i $hc_vcf -c haplotypecaller -o $hc_vcf\_HC.maf0\n\n";
    `$PYTHON/python $Bin/maf/vcf2maf0.py -i $hc_vcf -c haplotypecaller -o $hc_vcf\_HC.maf0`;
}

## Quality Filtering for haplotype caller 
print "$PYTHON/python $Bin/maf/pA_qSomHC.py < $hc_vcf\_HC.maf0 > $hc_vcf\_HC.maf1\n\n";
`$PYTHON/python $Bin/maf/pA_qSomHC.py < $hc_vcf\_HC.maf0 > $hc_vcf\_HC.maf1`;

## Convert to TCGA MAF format
print "$PYTHON/python $Bin/maf/oldMAF2tcgaMAF.py hg19 $hc_vcf\_HC.maf1 $hc_vcf\_HC.maf2\n\n";
`$PYTHON/python $Bin/maf/oldMAF2tcgaMAF.py hg19 $hc_vcf\_HC.maf1 $hc_vcf\_HC.maf2`;

## Select for Indels Only
print "$PYTHON/python $Bin/maf/indelOnly.py < $hc_vcf\_HC.maf2 > $hc_vcf\_qSomHC_InDels_TCGA_MAF.txt\n\n";
`$PYTHON/python $Bin/maf/indelOnly.py < $hc_vcf\_HC.maf2 > $hc_vcf\_qSomHC_InDels_TCGA_MAF.txt`;



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
    print "$PYTHON/python $Bin/maf/vcf2maf0.py -i $linked_vcf -aF $linked_txt -c mutect -o $linked_vcf\_maf0.txt -p $pairing -n $normal -t $tumor\n\n";
    `$PYTHON/python $Bin/maf/vcf2maf0.py -i $linked_vcf -aF $linked_txt -c mutect -o $linked_vcf\_maf0.txt -p $pairing -n $normal -t $tumor`;
    
    ## rescue 
    print "$PYTHON/python $Bin/rescue/DMP_rescue.py  < $linked_vcf\_maf0.txt > $linked_vcf\_maf1.txt 2> $linked_vcf\_maf_rescue.log\n\n";
    `$PYTHON/python $Bin/rescue/DMP_rescue.py  < $linked_vcf\_maf0.txt > $linked_vcf\_maf1.txt 2> $linked_vcf\_maf_rescue.log `;
    
    ## Change to tcga maf
    print "$PYTHON/python $Bin/maf/oldMAF2tcgaMAF.py hg19 $linked_vcf\_maf1.txt $linked_vcf\_maf2.txt\n\n";
    `$PYTHON/python $Bin/maf/oldMAF2tcgaMAF.py hg19 $linked_vcf\_maf1.txt $linked_vcf\_maf2.txt`;
    
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
`cut -f-39 $hc_vcf\_qSomHC_InDels_TCGA_MAF.txt > $output/$pre\_merge_maf0.txt`;

## For each mutect result, get all results except for the header
for my $vcf (@mutect_vcfs){
    my $linked_vcf = "$output/$vcf";
    print "tail -n +2 $linked_vcf\_DMPFilter_TCGA_MAF.txt >> $output/$pre\_merge_maf0.txt\n\n";
    `tail -n +2 $linked_vcf\_DMPFilter_TCGA_MAF.txt >> $output/$pre\_merge_maf0.txt`;
}

## link oncotator and oncotate!
symlink("$ONCOTATOR/db.properties", "./db.properties");
print "$ONCOTATOR/oncotateMaf.sh $output/$pre\_merge_maf0.txt $output/$pre\_merge_maf1.txt\n\n";
`$ONCOTATOR/oncotateMaf.sh $output/$pre\_merge_maf0.txt $output/$pre\_merge_maf1.txt`;

## fix Hugo?
print "$PYTHON/python $Bin/maf/pA_fixHugo.py <$output/$pre\_merge_maf1.txt >$output/$pre\_merge_maf2.txt\n\n";
`$PYTHON/python $Bin/maf/pA_fixHugo.py <$output/$pre\_merge_maf1.txt >$output/$pre\_merge_maf2.txt`;

## Update hugo symbol
print "$PERL/perl $Bin/update_gene_names_and_ids.pl $output/$pre\_merge_maf2.txt\n\n";
#`/bin/rm $Bin/lib/hugo_data.tsv`;
`$PERL/perl $Bin/update_gene_names_and_ids.pl $output/$pre\_merge_maf2.txt > $output/$pre\_merge_maf2_hugo.log 2>&1`;

## Grep for functional events YO!
print "cat $output/$pre\_merge_maf2.txt_hugo_modified | $PYTHON/python $Bin/maf/pA_Functional_Oncotator2.py > $output/$pre\_haplotect_TCGA_MAF.txt\n\n";
`cat $output/$pre\_merge_maf2.txt_hugo_modified | $PYTHON/python $Bin/maf/pA_Functional_Oncotator2.py > $output/$pre\_haplotect_TCGA_MAF.txt`;

## creating MAF for cbio portal submission
print "cat $output/$pre\_haplotect_TCGA_MAF.txt | $PYTHON/python $Bin/maf/pA_reSortCols.py $Bin/maf/finalCols_PORTAL.txt > $output/$pre\_haplotect_PORTAL_MAF.txt\n\n";
`cat $output/$pre\_haplotect_TCGA_MAF.txt | $PYTHON/python $Bin/maf/pA_reSortCols.py $Bin/maf/finalCols_PORTAL.txt > $output/$pre\_haplotect_PORTAL_MAF.txt`;

###       NOTE: This is different than the MAF4.txt because there are different headers here.
###       header information is missing everything after the first 39 columns, which is caller specific
print "creating clean MAF file\n";
print "/bin/cat $output/$pre\_merge_maf2.txt_hugo_modified | $PYTHON/python $Bin/maf/pA_Functional_Oncotator2.py | $PYTHON/python $Bin/maf/pA_reSortCols.py $Bin/maf/haplotect_abbrev_cols.txt > $output/$pre\_haplotect_merged_MAF.txt\n\n";
# Create nice MAF with essential columns
`/bin/cat $output/$pre\_merge_maf2.txt_hugo_modified | $PYTHON/python $Bin/maf/pA_Functional_Oncotator2.py | $PYTHON/python $Bin/maf/pA_reSortCols.py $Bin/maf/haplotect_abbrev_cols.txt > $output/$pre\_haplotect_merged_MAF.txt`;

print "annotating with cosmic\n";
print "$PYTHON/python $Bin/maf/maf_annotations/addCosmicAnnotation.py -i $output/$pre\_haplotect_TCGA_MAF.txt -o $output/$pre\_haplotect_TCGA_MAF_COSMIC_STANDARD.txt -f $Bin/data/CosmicMutantExport_v67_241013.tsv\n\n";
print "$PYTHON/python $Bin/maf/maf_annotations/addCosmicAnnotation.py -i $output/$pre\_haplotect_TCGA_MAF.txt -o $output/$pre\_haplotect_TCGA_MAF_COSMIC_DETAILED.txt -f $Bin/data/CosmicMutantExport_v67_241013.tsv -d\n\n";
`$PYTHON/python $Bin/maf/maf_annotations/addCosmicAnnotation.py -i $output/$pre\_haplotect_TCGA_MAF.txt -o $output/$pre\_haplotect_TCGA_MAF_COSMIC_STANDARD.txt -f $Bin/data/CosmicMutantExport_v67_241013.tsv`;
`$PYTHON/python $Bin/maf/maf_annotations/addCosmicAnnotation.py -i $output/$pre\_haplotect_TCGA_MAF.txt -o $output/$pre\_haplotect_TCGA_MAF_COSMIC_DETAILED.txt -f $Bin/data/CosmicMutantExport_v67_241013.tsv -d`;

print "annotating with mutation assessor\n";
print "$PYTHON/python $Bin/maf/maf_annotations/addMAannotation.py -i $output/$pre\_haplotect_TCGA_MAF_COSMIC_STANDARD.txt -o $output/$pre\_haplotect_TCGA_MAF_COSMIC_MA_STANDARD.txt\n";
print "$PYTHON/python $Bin/maf/maf_annotations/addMAannotation.py -i $output/$pre\_haplotect_TCGA_MAF_COSMIC_DETAILED.txt -o $output/$pre\_haplotect_TCGA_MAF_COSMIC_MA_DETAILED.txt -d\n\n";
`$PYTHON/python $Bin/maf/maf_annotations/addMAannotation.py -i $output/$pre\_haplotect_TCGA_MAF_COSMIC_STANDARD.txt -o $output/$pre\_haplotect_TCGA_MAF_COSMIC_MA_STANDARD.txt`;
`$PYTHON/python $Bin/maf/maf_annotations/addMAannotation.py -i $output/$pre\_haplotect_TCGA_MAF_COSMIC_DETAILED.txt -o $output/$pre\_haplotect_TCGA_MAF_COSMIC_MA_DETAILED.txt -d`;

if($patient && $bam_dir){
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

    print "$Bin/maf/fillout/GetBaseCountsMutliSample/GetBaseCountsMultiSample --fasta $REF_FASTA $bam_inputs --output $output/$pre/_haplotect_TCGA_basecounts.txt --maf $output/$pre\_haplotect_TCGA_MAF.txt --filter_improper_pair 0\n\n";
    `$Bin/maf/fillout/GetBaseCountsMutliSample/GetBaseCountsMultiSample --fasta $REF_FASTA $bam_inputs --output $output/$pre\_haplotect_TCGA_basecounts.txt --maf $output/$pre\_haplotect_TCGA_MAF.txt --filter_improper_pair 0`;

    print "$PYTHON/python $Bin/maf/fillout/dmp2portalMAF -m $output/$pre\_haplotect_TCGA_MAF.txt -p $pairing -P $patient -c haplotect -b $output/$pre\_haplotect_TCGA_basecounts.txt -o $output/$pre\_haplotect_TCGA_MAF_fillout.txt\n";
    `$PYTHON/python $Bin/maf/fillout/dmp2portalMAF -m $output/$pre\_haplotect_TCGA_MAF.txt -p $pairing -P $patient -c haplotect -b $output/$pre\_haplotect_TCGA_basecounts.txt -o $output/$pre\_haplotect_TCGA_MAF_fillout.txt`;
}


##
##
## Delete Temp Files
##
if($delete_temp){
    my @files = glob( $output . '/*.maf[012]');
    for my $fname (@files) {
        print "Removing: $fname \n";
        unlink($fname);
    }
    @files = glob ("$output/*mutect_calls* $output/*Haplotype* $output/*maf2.txt*");
    for my $fname (@files) {
        print "Removing: $fname \n";
        unlink($fname);
    }
}



