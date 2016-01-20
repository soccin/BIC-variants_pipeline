#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin);
use File::Path qw(make_path remove_tree);
use File::Basename;

my ($input, $somatic, $fastq, $target, $species, $output, $config, $help, $targetName);

GetOptions ('in_maf=s' => \$input,
'species=s' => \$species,
'output=s' => \$output,
'config=s' => \$config,
'somatic=s' => \$somatic,
'fastq|f' => \$fastq,
'target=s' => \$target,
'targetname=s' => \$targetName,
'help|h' => \$help ) or exit(1);

my $uID = `/usr/bin/id -u -n`;
chomp $uID;

if(!$input || !$species || !$output || !$config ){
    die "you are missing some input";
}

if($somatic){
$somatic = "_$somatic";
}else{
$somatic = "";
}

if($target){
    if(! -e $target){
        die "Target file $target does not exist: $!";
    }
    if(! $targetName){
        die "Target Name must be filled in if using -target option\n";
    }
}

my $PERL = '';
my $HG19_FASTA = '';
my $MM10_FASTA = '';
my $B37_FASTA = '';
my $BEDTOOLS = '';

open(CONFIG, "$config") or warn "CAN'T OPEN CONFIG FILE $config SO USING DEFAULT SETTINGS";
while(<CONFIG>){
    chomp;
    
    my @conf = split(/\s+/, $_);
    if($conf[0] =~ /perl/i){
        if(!-e "$conf[1]/perl"){
            die "CAN'T FIND perl IN $conf[1] $!";
        }
        $PERL = $conf[1];
    }
    elsif($conf[0] =~/bedtools/i){
        if(!-e "$conf[1]/bedtools"){
            die "CAN'T FIND bedtools in $conf[1] $!";
        }
        $BEDTOOLS = $conf[1];
    }
    elsif($conf[0] =~ /samtools/i){
        if(!-e "$conf[1]/samtools"){
            die "CAN'T FIND samtools IN $conf[1] $!";
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
    elsif($conf[0] =~ /mm10_fasta/i){
        if(!-e "$conf[1]"){
            die "CAN'T FIND $conf[1] $!";
        }
        $MM10_FASTA = $conf[1];
    }
    elsif($conf[0] =~ /b37_fasta/i){
        if(!-e "$conf[1]"){
            if($species =~ /human|^b37$/i){
                die "CAN'T FIND $conf[1] $!";
            }
        }
        $B37_FASTA = $conf[1];
    }
}
close CONFIG;

## FOR RIGHT NOW, ONLY HG19 IS BEING USED

my $REF_FASTA = $HG19_FASTA;
if($species =~ /b37/i){
    $REF_FASTA = $B37_FASTA;
}
if($species =~ /mm10/i){
    $REF_FASTA = $MM10_FASTA;
}
my $REF_DICT= $REF_FASTA;
$REF_DICT =~ s/\.[^.]+$//;
$REF_DICT = "$REF_DICT.dict";

# check to make sure config is not empty!
if(!-e $PERL) {
    die "Need PERL specified in config. $!";
}
if(!-e $BEDTOOLS){
    die "Need bedtools in config file. $!";
}

my $ref_base = basename($REF_FASTA);
# softlink reference
if(!-e "$output/ref$somatic"){
    mkdir("$output/ref$somatic", 0755) or die "Making ref didn't work $!";
}
if(!-e "$output/ref$somatic/$ref_base"){
    symlink($REF_FASTA, "$output/ref$somatic/$ref_base");
    symlink("$REF_FASTA.fai", "$output/ref$somatic/$ref_base.fai");
}

my $refDictBase = basename($REF_DICT);
if (!-e "$output/ref$somatic/$refDictBase"){
    symlink($REF_DICT, "$output/ref$somatic/$refDictBase");
}

if( ! -d "$output/bed$somatic/" ){
    print "$output/bed$somatic/ does not exist. Will create it now\n";
    mkdir("$output/bed$somatic", 0755) or die "Making tmp didn't work $!";
}
else{
    remove_tree("$output/bed$somatic") or die "Not able to remove directory $!";
    mkdir("$output/bed$somatic", 0755) or die "Making tmp didn't work $!";
}

# first no matter what make the maf a simple bed
open(my $in_fh, "<", $input);
open(my $out_fh, ">", "$output/bed$somatic/simple_maf.bed");
my @lines = <$in_fh>;
foreach my $line (@lines){
    #print "$line\n";
    next if ($line =~ m/^#|^Hugo_Symbol/);
    my @splitLine = split("\t", $line);
    my $zeroBase = $splitLine[5] - 1 ;
    print $out_fh "$splitLine[4]\t$zeroBase\t$splitLine[6]\n";
}
close $in_fh;
close $out_fh;


# pick fastq
if($fastq){
# take the hg19 Dict, and parse it into a genome file that bedtools takes in
    `grep "^\@SQ" $output/ref$somatic/$refDictBase | cut -f 2-3 | sed -e 's/SN://g' -e 's/LN://g' > $output/bed$somatic/bedtools_genome.txt`; 
    `$BEDTOOLS/bedtools slop -g $output/bed$somatic/bedtools_genome.txt -b 1 -i $output/bed$somatic/simple_maf.bed | $BEDTOOLS/bedtools getfasta -tab -fi $output/ref$somatic/$ref_base -fo $output/bed$somatic/triNucleotide.seq -bed -`;

    `mv $output/bed$somatic/triNucleotide.seq $output/triNucleotide.seq`;
}

if($target){
    `$BEDTOOLS/bedtools intersect -a $output/bed$somatic/simple_maf.bed -b $target -wa | $BEDTOOLS/bedtools sort -i - | uniq > $output/bed$somatic/temp_intersect.txt`;

    open(my $in_temp, "<", "$output/bed$somatic/temp_intersect.txt");
    open(my $out_temp, ">", "$output/bed$somatic/maf_targets.$targetName"); 
    my @lines = <$in_temp>;
    foreach my $line (@lines){
        my @parts = split("\t", $line);
        my $start = $parts[1] + 1;
        print $out_temp "$parts[0]:$start-$parts[2]\n";
    }

    close $in_temp;
    close $out_temp;

    `mv $output/bed$somatic/maf_targets.$targetName $output/maf_targets.$targetName`;
}

