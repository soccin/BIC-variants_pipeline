#!/opt/bin/perl

use strict;
use Getopt::Long qw(GetOptions);

my ($pair, $vcf, $pre, $Bin, $help, $outdir);
my $uID = `/usr/bin/id -u -n`;
chomp $uID;


$pre = 'TEMP';
GetOptions ('pair=s' => \$pair,
            'vcf=s' => \$vcf,
            'pre=s' => \$pre,
            'bin=s' => \$Bin,
            'outdir=s' => \$outdir,
            'help' => \$help) or exit(1);

if(!$pair || $help){
    print <<HELP;

    USAGE: ./pairing_qc.pl -vcf VCF -pair PAIR -pre PRE 
        * PAIR: file listing tumor/normal pairing of samples for mutect/maf conversion; if not specified, considered unpaired
        * VCF: population VCF file from which SNPs will be selected for pair evaluation
        * PRE: output prefix (default: TEMP)
        * BIN: path to root of exome pipeline
        * OUTDIR: output directory for all fingerprint files
HELP
exit;
}

my $alerts_dir = "/ifs/data/alerts/$pre";

## move to output dir
`cd $outdir`;

## make a copy of pairing file in out directory for web app
`cp $pair $outdir`;

## filter VCF
`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_FP -hold_jid $pre\_$uID\_HC -q lau.q $Bin/qCMD $Bin/qc/fingerprint.py $vcf $outdir/$pre\_FINGERPRINT.vcf`;

## flatten & further filter VCF
`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_FFP -hold_jid $pre\_$uID\_FP -q lau.q $Bin/qCMD "cat $outdir/$pre\_FINGERPRINT.vcf | python $Bin/qc/flattenFingerprint.py > $outdir/$pre\_FINGERPRINT_FLAT.txt"`;
## plot flattened fingerprint
`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_PLOT_FFP -hold_jid $pre\_$uID\_FFP -q lau.q $Bin/qCMD /opt/common/R/R-3.0.3/bin/Rscript $Bin/qc/plotFingerprint.R $outdir/$pre\_FINGERPRINT_FLAT.txt`;


## pairing by pval
`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_FINDPAIRS_PVAL -hold_jid $pre\_$uID\_FP -q lau.q $Bin/qCMD $Bin/qc/pval_findPairs.py -vcf $outdir/$pre\_FINGERPRINT.vcf -out $outdir/$pre\_FINGERPRINT.vcf_pval_pairs.txt`;
### clustering by pval
`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_CLUSTER_PVAL -hold_jid $pre\_$uID\_FINDPAIRS_PVAL -q lau.q $Bin/qCMD $Bin/qc/pval_cluster.py $outdir/$pre\_FINGERPRINT.vcf_pval_pairs.txt $outdir/$pre\_FINGERPRINT.vcf_pval_pairs_cluster.txt 0.05`;

### pairing by IBS
`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_FINDPAIRS_IBS -hold_jid $pre\_$uID\_FP -q lau.q $Bin/qCMD $Bin/qc/ibs_findPairs.py $outdir/$pre\_FINGERPRINT.vcf $outdir/$pre\_FINGERPRINT.vcf_ibs_pairs.txt`;
### clustering by IBS
`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_CLUSTER_IBS -hold_jid $pre\_$uID\_FINDPAIRS_IBS -q lau.q $Bin/qCMD $Bin/qc/ibs_cluster.py $outdir/$pre\_FINGERPRINT.vcf_ibs_pairs.txt $outdir/$pre\_FINGERPRINT.vcf_ibs_pairs_cluster.txt 1.8 0.3`;
### IBS plot
`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_PLOT_IBS -hold_jid $pre\_$uID\_FINDPAIRS_IBS -q lau.q $Bin/qCMD /opt/common/R/R-3.0.3/bin/Rscript $Bin/qc/ibs_plot.R \"\\\"projID='$pre'\\\"\" \"\\\"manifestPairing='$pair'\\\"\" \"\\\"IBSFile='$outdir/$pre\_FINGERPRINT.vcf_ibs_pairs.txt'\\\"\"`;


#`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_ALERT -hold_jid $pre\_$uID\_CLUSTER_IBS,$pre\_$uID\_CLUSTER_PVAL,$pre\_$uID\_FFP -q lau.q $Bin/qCMD $Bin/qc/alert.sh $pair $outdir $alerts_dir`;
