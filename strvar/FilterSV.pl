#!/usr/bin/perl -w
##########FilterSV.pl########
#Author: Ronak Shah
#Date: 01/02/2014
#Version: 1.0
#Description: Filter the VCF file for somatic SV's
#process.
#############
###Log:
##01/02/2014
#v1.0
use lib qw(/ifs/data/zeng/dmp/resources/perllib);
use strict;
use Getopt::Long;
use POSIX;

#use Statistics::Test::WilcoxonRankSum;
#use Statistics::DependantTTest;
#use Statistics::Distributions;
my (
	$vcf,                              $hotspotFile,
	$tumorId,                          $normalId,
	$outdir,                           $outputFile,
	$overallSupportingReadsPE,         $sampleTumorSupportingReadsPE,
	$sampleNormalSupportingReadsPE,    $overallSupportingReadsSR,
	$overallSupportingReadsPE_Hotspot, $sampleTumorSupportingReadsPE_Hotspot,
	$sampleNormalSupportingReadsPE_Hotspot, $overallSupportingReadsSR_Hotspot,
	$lengthOfSV,                            $overallMAPQ,
	$overallMAPQ_Hotspot,                   $sampleTumorGenotypeQualityFilter,
	$sampleTumorGenotypeQualityFilterHotspot,
	$sampleNormalGenotypeQualityFilter,
	$sampleNormalGenotypeQualityFilterHotspot,    #$rcProbabilityHotspot,
	$tnRatioHotspot,                              #$rcProbability,
	$tnRatio, $annotationInput, $TypeOfSV
);

#--This variable holds the current time
my $now = time;
if (
	@ARGV < 1
	or !GetOptions(
		'vcf|i:s'                        => \$vcf,
		'hotspotFile|hsf:s'              => \$hotspotFile,
		'typeOfSV|tos:s'                 => \$TypeOfSV,
		'tumorId|tid:s'                  => \$tumorId,
		'normalId|nid:s'                 => \$normalId,
		'outdir:s'                       => \$outdir,
		'outputFile|o:s'                 => \$outputFile,
		'overallSupportingReadsPE|ope:i' => \$overallSupportingReadsPE,
		'overallSupportingReadsSR|osr:i' => \$overallSupportingReadsSR,
		'overallSupportingReadsPEforHotspot|opeh:i' =>
		  \$overallSupportingReadsPE_Hotspot,
		'overallSupportingReadsSRforHotspot|osrh:i' =>
		  \$overallSupportingReadsSR_Hotspot,
		'sampleTumorSupportingReadsPE|stpe:i' => \$sampleTumorSupportingReadsPE,
		'sampleNormalSupportingReadsPE|snpe:i' =>
		  \$sampleNormalSupportingReadsPE,
		'sampleTumorSupportingReadsPEforHotspot|stpeh:i' =>
		  \$sampleTumorSupportingReadsPE_Hotspot,
		'sampleNormalSupportingReadsPEforHotspot|snpeh:i' =>
		  \$sampleNormalSupportingReadsPE_Hotspot,
		'lengthOfSV|svlen:i'           => \$lengthOfSV,
		'overallMAPQ|omq:i'            => \$overallMAPQ,
		'overallMAPQForHotspot|omqh:i' => \$overallMAPQ_Hotspot,
		'sampleTumorGenotypeQualityFilter|stgqf:i' =>
		  \$sampleTumorGenotypeQualityFilter,
		'sampleNormalGenotypeQualityFilter|sngqf:i' =>
		  \$sampleNormalGenotypeQualityFilter,
		'sampleTumorGenotypeQualityFilterHotspot|stgqfh:i' =>
		  \$sampleTumorGenotypeQualityFilterHotspot,
		'sampleNormalGenotypeQualityFilterHotspot|sngqfh:i' =>
		  \$sampleNormalGenotypeQualityFilterHotspot,
		'TNratio|tnr:i'         => \$tnRatio,
		'TNratioHotspot|tnrh:i' => \$tnRatioHotspot

		  #'RCprobability|rcr:i'        => \$rcProbability,
		  #'RCpobabilityHotspot|rcrh:i' => \$rcProbabilityHotspot
	)
  )
{
	Usage();
}

#Checks for Input
#vcf
if ( ( !$vcf ) or ( !( -f $vcf ) ) ) {
	print "Please enter the name of the vcf file to be analyzed. See Usage.\n";
	Usage();
	exit(1);
}
else {
	print "VCF:$vcf\n";
}

#hotspotFile
if ( ( !$hotspotFile ) or ( !( -f $hotspotFile ) ) ) {
	print
	  "Please enter the name of the hotspot file to be analyzed. See Usage.\n";
	Usage();
	exit(1);
}
else {
	print "HotspotFile:$hotspotFile\n";
}

#tumorId
if ( !$tumorId ) {
	print "Please enter the id for the tumor bam file. See Usage.\n";
	Usage();
	exit(1);
}
else {
	print "TumorID:$tumorId\n";
}

#normalId
if ( !$normalId ) {
	print "Please enter the id for the normal bam file. See Usage.\n";
	Usage();
	exit(1);
}
else {
	print "NormalID:$normalId\n";
}

#TypeOfSV
if ( !$TypeOfSV ) {
	print "Please enter the type of Structural Variant. See Usage.\n";
	Usage();
	exit(1);
}
else {
	print "TypeOfSV:$TypeOfSV\n";
}

#outdir
if ( !$outdir ) {
	print "Output dir not given default will be used.\n";
	$outdir = getcwd;
	print "OUPUTDIR:$outdir\n";
}
else {
	if ( !( -d $outdir ) ) {
		print "Output dir given is incorrect or doesnot exists.\n";
		exit(1);
	}
	else {
		print "OUTPUTDIR:$outdir.\n";
	}
}

#outputFile
if ( !$outputFile ) {
	print "Please enter the output file name. See Usage.\n";
	Usage();
	exit(1);
}
else {
	if ( $outputFile =~ /\./ ) {
		($annotationInput) = $outputFile =~ /(.*)\./;
		$annotationInput = $annotationInput . "_dRangerInput.txt";
	}
	else {
		($annotationInput) = $outputFile;
		$annotationInput = $annotationInput . "_dRangerInput.txt";
	}
}

#overallSupportingReadsPE
if ( !$overallSupportingReadsPE ) {
	print "overallSupportingReadsPE not given default will be used.\n";
	$overallSupportingReadsPE = 5;
	print "overallSupportingReadsPE:$overallSupportingReadsPE\n";
}
else {
	print "overallSupportingReadsPE:$overallSupportingReadsPE.\n";
}

#sampleTumorSupportingReadsPE
if ( !$sampleTumorSupportingReadsPE ) {
	print "sampleTumorSupportingReadsPE not given default will be used.\n";
	$sampleTumorSupportingReadsPE = 2;
	print "sampleTumorSupportingReadsPE:$sampleTumorSupportingReadsPE\n";
}
else {
	print "sampleTumorSupportingReadsPE:$sampleTumorSupportingReadsPE.\n";
}

#sampleNormalSupportingReadsPE
if ( !$sampleNormalSupportingReadsPE ) {
	print "sampleNormalSupportingReadsPE not given default will be used.\n";
	$sampleNormalSupportingReadsPE = 2;
	print "sampleNormalSupportingReadsPE:$sampleNormalSupportingReadsPE\n";
}
else {
	print "sampleNormalSupportingReadsPE:$sampleNormalSupportingReadsPE.\n";
}

#sampleSupportingReadsPE_Hotspot
if ( !$sampleTumorSupportingReadsPE_Hotspot ) {
	print
	  "sampleTumorSupportingReadsPE_Hotspot not given default will be used.\n";
	$sampleTumorSupportingReadsPE_Hotspot = 1;
	print
"sampleTumorSupportingReadsPE_Hotspot:$sampleTumorSupportingReadsPE_Hotspot\n";
}
else {
	print
"sampleTumorSupportingReadsPE_Hotspot:$sampleTumorSupportingReadsPE_Hotspot.\n";
}

#sampleNormalSupportingReadsPE_Hotspot
if ( !$sampleNormalSupportingReadsPE_Hotspot ) {
	print
	  "sampleNormalSupportingReadsPE_Hotspot not given default will be used.\n";
	$sampleNormalSupportingReadsPE_Hotspot = 3;
	print
"sampleNormalSupportingReadsPE_Hotspot:$sampleNormalSupportingReadsPE_Hotspot\n";
}
else {
	print
"sampleNormalSupportingReadsPE_Hotspot:$sampleNormalSupportingReadsPE_Hotspot.\n";
}

#overallSupportingReadsSR
if ( !$overallSupportingReadsSR ) {
	print "overallSupportingReadsSR not given default will be used.\n";
	$overallSupportingReadsSR = 0;
	print "overallSupportingReadsSR:$overallSupportingReadsSR\n";
}
else {
	print "overallSupportingReadsSR:$overallSupportingReadsSR.\n";
}

#overallSupportingReadsPE_Hotspot
if ( !$overallSupportingReadsPE_Hotspot ) {
	print "overallSupportingReadsPE_Hotspot not given default will be used.\n";
	$overallSupportingReadsPE_Hotspot = 3;
	print
	  "overallSupportingReadsPE_Hotspot:$overallSupportingReadsPE_Hotspot\n";
}
else {
	print
	  "overallSupportingReadsPE_Hotspot:$overallSupportingReadsPE_Hotspot.\n";
}

#overallSupportingReadsSR_Hotspot
if ( !$overallSupportingReadsSR_Hotspot ) {
	print "overallSupportingReadsSR_Hotspot not given default will be used.\n";
	$overallSupportingReadsSR_Hotspot = 0;
	print
	  "overallSupportingReadsSR_Hotspot:$overallSupportingReadsSR_Hotspot\n";
}
else {
	print
	  "overallSupportingReadsSR_Hotspot:$overallSupportingReadsSR_Hotspot.\n";
}

#lengthOfSV
if ( !$lengthOfSV ) {
	print "lengthOfSV not given default will be used.\n";
	$lengthOfSV = 500;
	print "lengthOfSV:$lengthOfSV\n";
}
else {
	print "lengthOfSV:$lengthOfSV.\n";
}

#overallMAPQ
if ( !$overallMAPQ ) {
	print "overallMAPQ not given default will be used.\n";
	$overallMAPQ = 10;
	print "overallMAPQ:$overallMAPQ\n";
}
else {
	print "overallMAPQ:$overallMAPQ.\n";
}

#overallMAPQ_Hotspot
if ( !$overallMAPQ_Hotspot ) {
	print "overallMAPQ_Hotspot not given default will be used.\n";
	$overallMAPQ_Hotspot = 5;
	print "overallMAPQ_Hotspot:$overallMAPQ_Hotspot\n";
}
else {
	print "overallMAPQ_Hotspot:$overallMAPQ_Hotspot.\n";
}

#sampleTumorGenotypeQualityFilter
if ( !$sampleTumorGenotypeQualityFilter ) {
	print "sampleTumorGenotypeQualityFilter not given default will be used.\n";
	$sampleTumorGenotypeQualityFilter = 15;
	print
	  "sampleTumorGenotypeQualityFilter = $sampleTumorGenotypeQualityFilter\n";
}
else {
	print
	  "sampleTumorGenotypeQualityFilter = $sampleTumorGenotypeQualityFilter\n";
}

#sampleNormalGenotypeQualityFilter
if ( !$sampleNormalGenotypeQualityFilter ) {
	print "sampleNormalGenotypeQualityFilter not given default will be used.\n";
	$sampleNormalGenotypeQualityFilter = 15;
	print
"sampleNormalGenotypeQualityFilter = $sampleNormalGenotypeQualityFilter\n";
}
else {
	print
"sampleNormalGenotypeQualityFilter = $sampleNormalGenotypeQualityFilter\n";
}

#sampleTumorGenotypeQualityFilterHotspot
if ( !$sampleTumorGenotypeQualityFilterHotspot ) {
	print
"sampleTumorGenotypeQualityFilterHotspot not given default will be used.\n";
	$sampleTumorGenotypeQualityFilterHotspot = 5;
	print
"sampleTumorGenotypeQualityFilterHotspot = $sampleTumorGenotypeQualityFilterHotspot\n";
}
else {
	print
"sampleTumorGenotypeQualityFilterHotspot = $sampleTumorGenotypeQualityFilterHotspot\n";
}

#sampleGenotypeQualityFilterHotspot
if ( !$sampleNormalGenotypeQualityFilterHotspot ) {
	print
"sampleNormalGenotypeQualityFilterHotspot not given default will be used.\n";
	$sampleNormalGenotypeQualityFilterHotspot = 20;
	print
"sampleNormalGenotypeQualityFilterHotspot = $sampleNormalGenotypeQualityFilterHotspot\n";
}
else {
	print
"sampleNormalGenotypeQualityFilterHotspot = $sampleNormalGenotypeQualityFilterHotspot\n";
}

#tnRatio
if ( !$tnRatio ) {
	print "tnRatio not given default will be used.\n";
	$tnRatio = 0;
	print "tnRatio:$tnRatio\n";
}
else {
	print "tnRatio:$tnRatio.\n";
}

#tnRatioHotspot
if ( !$tnRatioHotspot ) {
	print "tnRatioHotspot not given default will be used.\n";
	$tnRatioHotspot = 0;
	print "tnRatioHotspot:$tnRatioHotspot\n";
}
else {
	print "tnRatioHotspot:$tnRatioHotspot.\n";
}

=begin
#rcRatio
if ( !$rcProbability ) {
	print "rcProbability not given default will be used.\n";
	$rcProbability = 0.01;
	print "rcProbability:$rcProbability\n";
}
else {
	print "rcProbability:$rcProbability.\n";
}

#rcRatioHotspot
if ( !$rcProbabilityHotspot ) {
	print "rcProbabilityHotspot not given default will be used.\n";
	$rcProbabilityHotspot = 0.1;
	print "rcProbabilityHotspot:$rcProbabilityHotspot\n";
}
else {
	print "rcProbabilityHotspot:$rcProbabilityHotspot.\n";
}
=cut

#Run Filtering of VCF
my ( $tumorVals, $normalVals, $info ) = &FilterVcf();

#--Calculate total runtime (current time minus start time)
$now = time - $now;

#--Print runtime
printf(
	"\n\nTotal running time: %02d:%02d:%02d\n\n",
	int( $now / 3600 ),
	int( ( $now % 3600 ) / 60 ),
	int( $now % 60 )
);
print "\n!!!!Done, Thanks for using the script!!!\n";
exit;
#####################################
#####################################
#How to use the script.
sub Usage {
	print "Unknow option: @_\n" if (@_);
	print "\nUsage : FilterSV.pl [options]
        [--vcf|i 						S	 	1000 Genomes Format Structural Variant[SV] VCF.]
        [--hotspotFile|hsf 						S	 	BED4 Format bed file.]
        [--typeOfSV|tos 						S	 	Type of structural Variant(please use these terms accordingly: Deletion:DEL,Duplication:DUP,Inversion:INV,Translocation:TRA)]
        [--tumorId|tid					S		ID for the tumor bam.]
        [--normalId|nid					S		ID for the Noraml bam.]
        [--outdir						S		Output directory.]
        [--outputFile|o					S		Name of the Output file.]
        [--overallSupportingReadsPE|pe			I		Threshold for number of pair-end reads supporting the SV.]
        [--overallSupportingReadsSR|sr			I		Threshold for number of split reads supporting the SV.]
        [--lengthOfSV|svlen	I	Threshold for length of the SV.]
        [--overallMAPQ|mq	I	Threshold for Mapping Quality.]
        [--sampleGenotypeQualityFilter|gqf	I	Genotype Quality Filter for non-hotspot region. default(15:This will eliminate all SV's that have GQ < 15 in Tumor & GQ > 15 in Normal)]
		[--sampleGenotypeQualityFilterHotspot|gqfh	I	Genotype Quality Filter for hotspot region. default(5:This will eliminate all SV's that have GQ < 5 in Tumor)]
		
	\n";
	print "For any question please contact Ronak Shah (shahr2\@mskcc.org)\n";
	print "!~!~! Comments, flames and suggestion are welcome !~!~!\n";
	exit;
}
#####################################
#####################################
#Read the Hotspot Locations from Bed4 File
sub ReadHotspotFile {
	my %hotspotLocs = ();
	open( HFH, $hotspotFile )
	  || die "Cannot open file $hotspotFile, Error: $!\n";
	while (<HFH>) {
		chomp($_);
		my ( $chr, $start, $end, $gene ) = split( "\t", $_ );

		#print "$chr, $start, $end, $gene\n";
		if ( exists $hotspotLocs{$chr} ) {
			my $val = $hotspotLocs{$chr};
			$hotspotLocs{$chr} = $val . "#" . $start . ":" . $end . ":" . $gene;
		}
		else {
			$hotspotLocs{$chr} = $start . ":" . $end . ":" . $gene;
		}
	}
	close(HFH);
	return ( \%hotspotLocs );
}
#####################################
#####################################
#Check if location is hotspot or not
sub CheckHotspot {
	my ( $chr1, $pos1, $chr2, $pos2 ) = @_;

	#Read Hotspot Locations
	my $hotspotHash = &ReadHotspotFile();
	my %hotspotLocs = %$hotspotHash;
	my $hotspotFlag;
	if ( exists $hotspotLocs{$chr1} ) {
		my $hotspotLocations = $hotspotLocs{$chr1};
		if ( $hotspotLocations =~ m/#/ ) {
			my @multiLocs = split( "#", $hotspotLocations );
			foreach my $hsLocation (@multiLocs) {
				my ( $hsStart, $hsEnd, $hsGene ) =
				  split( ":", $hsLocation );
				if ( $pos1 >= $hsStart && $pos1 <= $hsEnd ) {
					$hotspotFlag = 2;
					last;
				}
				else {
					$hotspotFlag = 1;
				}
			}
		}
		else {
			my ( $hsStart, $hsEnd, $hsGene ) =
			  split( ":", $hotspotLocations );
			if ( $pos1 >= $hsStart && $pos1 <= $hsEnd ) {
				$hotspotFlag = 2;
			}
			else {
				$hotspotFlag = 1;
			}
		}
	}
	elsif ( exists $hotspotLocs{$chr2} ) {
		my $hotspotLocations = $hotspotLocs{$chr2};
		if ( $hotspotLocations =~ m/#/ ) {
			my @multiLocs = split( "#", $hotspotLocations );
			foreach my $hsLocation (@multiLocs) {
				my ( $hsStart, $hsEnd, $hsGene ) =
				  split( ":", $hsLocation );
				if ( $pos2 >= $hsStart && $pos2 <= $hsEnd ) {
					$hotspotFlag = 2;
					last;
				}
				else {
					$hotspotFlag = 1;
				}
			}
		}
		else {
			my ( $hsStart, $hsEnd, $hsGene ) =
			  split( ":", $hotspotLocations );
			if ( $pos2 >= $hsStart && $pos2 <= $hsEnd ) {
				$hotspotFlag = 2;
			}
			else {
				$hotspotFlag = 1;
			}
		}
	}
	else {
		$hotspotFlag = 1;
	}
	return ($hotspotFlag);
}
#####################################
#####################################
#First Tier Filter for VCF file.
sub FilterVcf {
	my @header = ();
	my (
		$tumorIndex, $normalIndex, $infoIndex,
		$chrIndex,   $startIndex,  $filterIndex
	);

	#Open the Output File
	open( OFH, ">", "$outdir/$outputFile" )
	  || die "Cannot open file $outdir/$outputFile, Error:$!\n";

	#Open the dranger input File
	open( DFH, ">", "$outdir/$annotationInput" )
	  || die "Cannot open file $outdir/$annotationInput, Error:$!\n";
	print DFH "chr1\tpos1\tstr1\tchr2\tpos2\tstr2\n";

	#Open the VCF filet to filter
	open( VFH, $vcf ) || die "Cannot open file $vcf, Error:$!\n";
	while ( my $line = <VFH> ) {
		chomp($line);

		#Get the Header
		if ( $line =~ m/^#/ ) {
			print OFH "$line\n";

			#Get what is tumor what is normal
			if ( $line =~ m/^#CHROM/ ) {
				@header = split( "\t", $line );
				for ( my $i = 0 ; $i < scalar(@header) ; $i++ ) {
					if ( $header[$i] =~ /$tumorId/ ) {
						$tumorIndex = $i;
					}
					if ( $header[$i] =~ /$normalId/ ) {
						$normalIndex = $i;
					}
					if ( $header[$i] =~ /INFO/ ) {
						$infoIndex = $i;
					}
					if ( $header[$i] =~ /CHROM/ ) {
						$chrIndex = $i;
					}
					if ( $header[$i] =~ /POS/ ) {
						$startIndex = $i;
					}
					if ( $header[$i] =~ /FILTER/ ) {
						$filterIndex = $i;
					}
				}
			}
		}
		else {

			#Get All the Lines
			my (@line) = split( "\t", $line );

			#Get Chromosome
			my $svChr = $line[$chrIndex];
			next if ( $svChr =~ /MT/ );

			#Get Start Position
			my $svStart = $line[$startIndex];

			#Get Filter String
			my $filterString = $line[$filterIndex];

			#Get Info column
			my @info = split( ";", $line[$infoIndex] );
			my %infoData = ();
			foreach my $infoTypes (@info) {
				if ( $infoTypes =~ m/=/ ) {
					my ( $key, $value ) = split( "=", $infoTypes );
					$infoData{$key} = $value;
				}
				else {
					$infoData{$infoTypes} = "";
				}
			}

			#get Tumor Genotype Information
			my ( $tGT, $tGL, $tGQ, $tFT, $tRC, $tDR, $tDV,$tRR, $tRV ) =
			  split( ":", $line[$tumorIndex] );

			#Get Normal Genotype Information
			my ( $nGT, $nGL, $nGQ, $nFT, $nRC, $nDR, $nDV,$nRR, $nRV ) =
			  split( ":", $line[$normalIndex] );
			my ( $svlen, $mapq, $peReads, $srReads, $svType, $svTNratio,
				$svChr2, $svEnd )
			  = 0.0;
			my ( $connectionType, $str1, $str2 ) = "";

			#Calculate Tumor Noraml SV read ratio if its 0
			#($svTNratio) = 5 * $nDV;
			$svlen          = $infoData{"SVLEN"}  if exists $infoData{"SVLEN"};
			$mapq           = $infoData{"MAPQ"}   if exists $infoData{"MAPQ"};
			$svType         = $infoData{"SVTYPE"} if exists $infoData{"SVTYPE"};
			$peReads        = $infoData{"PE"}     if exists $infoData{"PE"};
			$srReads        = $infoData{"SR"}     if exists $infoData{"SR"};
			$svEnd          = $infoData{"END"}    if exists $infoData{"END"};
			$connectionType = $infoData{"CT"}     if exists $infoData{"CT"};
			$svChr2         = $infoData{"CHR2"}   if exists $infoData{"CHR2"};
			if ( !$svChr2 ) { $svChr2 = $svChr }
			my $convertedSVchr;
			$convertedSVchr = $svChr;
			if($convertedSVchr =~ m/^chr/ )
			{
				$convertedSVchr = substr($convertedSVchr, 3);	
			}
			$convertedSVchr = "23" if ( $svChr =~ /X/ );
			$convertedSVchr = "24" if ( $svChr =~ /Y/ );
			my $convertedSVchr2;
			$convertedSVchr2 = $svChr2;
			if($convertedSVchr2 =~ m/^chr/)
			{
				$convertedSVchr2 = substr($convertedSVchr2, 3);
			}
			$convertedSVchr2 = "23" if ( $svChr2 =~ /X/ );
			$convertedSVchr2 = "24" if ( $svChr2 =~ /Y/ );

			#Get Connection Info
			my ( $startCT, $endCT ) = split( "to", $connectionType );
			$str1 = "0" if ( $startCT == 3 );
			$str2 = "0" if ( $endCT == 3 );
			$str1 = "1" if ( $startCT == 5 );
			$str2 = "1" if ( $endCT == 5 );
			my $filterFlag;
			if ( $TypeOfSV =~ /DEL/i ) {
				$filterFlag = &FilterDeletion(
					$svChr,        $svStart, $svChr2, $svEnd,
					$filterString, $svlen,   $mapq,   $peReads,
					$srReads,      $tFT,     $nFT,    $tGQ,
					$nGQ,          $tDV,     $nDV
				);
			}
			if ( $TypeOfSV =~ /DUP/i ) {
				$filterFlag = &FilterDuplication(
					$svChr,        $svStart, $svChr2, $svEnd,
					$filterString, $svlen,   $mapq,   $peReads,
					$srReads,      $tFT,     $nFT,    $tGQ,
					$nGQ,          $tDV,     $nDV
				);
			}
			if ( $TypeOfSV =~ /INV/i ) {
				$filterFlag = &FilterInversion(
					$svChr,        $svStart, $svChr2, $svEnd,
					$filterString, $svlen,   $mapq,   $peReads,
					$srReads,      $tFT,     $nFT,    $tGQ,
					$nGQ,          $tDV,     $nDV
				);
			}
			if ( $TypeOfSV =~ /TRA/i ) {
				$filterFlag = &FilterTranslocation(
					$svChr,        $svStart, $svChr2, $svEnd,
					$filterString, $svlen,   $mapq,   $peReads,
					$srReads,      $tFT,     $nFT,    $tGQ,
					$nGQ,          $tDV,     $nDV
				);
			}

			#print "$hotspotFlag:$svChr\t$svStart\n";
			#Check Hotspot and filter
			if ( $filterFlag == 2 ) {
				print OFH "$line\n";
				print DFH
"$convertedSVchr\t$svStart\t$str1\t$convertedSVchr2\t$svEnd\t$str2\n";
			}
			else {
				next;
			}
		}
	}
	close(VFH);
	close(OFH);
	close(DFH);
	return;
}
#####################################
#####################################
#Filter record for deletions.
sub FilterDeletion {
	my (
		$svChr, $svStart, $svChr2,  $svEnd,   $filterString,
		$svlen, $mapq,    $peReads, $srReads, $tFT,
		$nFT,   $tGQ,     $nGQ,     $tDV,     $nDV
	) = @_;
	my $filterFlag = 1;

	#Start Checking if its Hotspot or not.
	my $hotspotFlag = CheckHotspot( $svChr, $svStart, $svChr2, $svEnd );
	if ( $hotspotFlag == 2 ) {

		#$svTNratio = $tnRatioHotspot if (! $svTNratio);
		my $dupFlag = 1;
		if (    ( $filterString =~ /PASS/ )
			and ( $nFT =~ /LowQual/ )
			and ( $svlen >= $lengthOfSV ) )
		{
			$filterFlag = 2;
			$dupFlag    = 2;
		}
		if (
			    ( $svlen >= $lengthOfSV )
			and ( $mapq >= $overallMAPQ_Hotspot )
			and ( $peReads >= $overallSupportingReadsPE_Hotspot )
			and ( $srReads >= $overallSupportingReadsSR_Hotspot )

			and ( $tDV >= $sampleTumorSupportingReadsPE_Hotspot )
			and ( $nDV <= $sampleNormalSupportingReadsPE_Hotspot )

			#and ( $svTNratio >= $tnRatioHotspot )
			#and ( $tGQ >= $sampleTumorGenotypeQualityFilterHotspot )
			#and ( $nGQ <= $sampleNormalGenotypeQualityFilter )
			and ( $dupFlag == 1 )
		  )
		{
			$filterFlag = 2;
		}
	}
	else {
		if (
			    ( $svlen >= $lengthOfSV )
			and ( $mapq >= $overallMAPQ )
			and ( $peReads >= $overallSupportingReadsPE )
			and ( $srReads >= $overallSupportingReadsSR )

			and ( $tDV >= $sampleTumorSupportingReadsPE )
			and ( $nDV <= $sampleNormalSupportingReadsPE )

			#and ( $svTNratio >= $tnRatio )
			#and ( $tGQ >= $sampleTumorGenotypeQualityFilter )
			#and ( $nGQ <= $sampleNormalGenotypeQualityFilter )
		  )
		{
			$filterFlag = 2;
		}
	}
	return ($filterFlag);
}
#####################################
#####################################
#Filter record for duplications
sub FilterDuplication {
	my (
		$svChr, $svStart, $svChr2,  $svEnd,   $filterString,
		$svlen, $mapq,    $peReads, $srReads, $tFT,
		$nFT,   $tGQ,     $nGQ,     $tDV,     $nDV
	) = @_;
	my $filterFlag = 1;

	#Start Checking if its Hotspot or not.
	my $hotspotFlag = CheckHotspot( $svChr, $svStart, $svChr2, $svEnd );
	if ( $hotspotFlag == 2 ) {

		#$svTNratio = $tnRatioHotspot if (! $svTNratio);
		my $dupFlag = 1;
		if (    ( $filterString =~ /PASS/ )
			and ( $nFT =~ /LowQual/ )
			and ( $svlen >= $lengthOfSV ) )
		{
			$filterFlag = 2;
			$dupFlag    = 2;
		}
		if (
			    ( $svlen >= $lengthOfSV )
			and ( $mapq >= $overallMAPQ_Hotspot )
			and ( $peReads >= $overallSupportingReadsPE_Hotspot )
			and ( $srReads >= $overallSupportingReadsSR_Hotspot )

			and ( $tDV >= $sampleTumorSupportingReadsPE_Hotspot )
			and ( $nDV <= $sampleNormalSupportingReadsPE_Hotspot )

			#and ( $svTNratio >= $tnRatioHotspot )
			#and ( $tGQ >= $sampleTumorGenotypeQualityFilterHotspot )
			#and ( $nGQ <= $sampleNormalGenotypeQualityFilter )
			and ( $dupFlag == 1 )
		  )
		{
			$filterFlag = 2;
		}
	}
	else {
		if (
			    ( $svlen >= $lengthOfSV )
			and ( $mapq >= $overallMAPQ )
			and ( $peReads >= $overallSupportingReadsPE )
			and ( $srReads >= $overallSupportingReadsSR )

			#and ( $svTNratio >= $tnRatio )
			and ( $tDV >= $sampleTumorSupportingReadsPE )
			and ( $nDV <= $sampleNormalSupportingReadsPE )
			#and ( $tGQ >= $sampleTumorGenotypeQualityFilter )
			#and ( $nGQ <= $sampleNormalGenotypeQualityFilter )
		  )
		{
			$filterFlag = 2;
		}
	}
	return ($filterFlag);
}
#####################################
#####################################
#Filter record for inversions
sub FilterInversion {
	my (
		$svChr, $svStart, $svChr2,  $svEnd,   $filterString,
		$svlen, $mapq,    $peReads, $srReads, $tFT,
		$nFT,   $tGQ,     $nGQ,     $tDV,     $nDV
	) = @_;
	my $filterFlag = 1;

	#Start Checking if its Hotspot or not.
	my $hotspotFlag = CheckHotspot( $svChr, $svStart, $svChr2, $svEnd );
	if ( $hotspotFlag == 2 ) {

		#$svTNratio = $tnRatioHotspot if (! $svTNratio);
		my $dupFlag = 1;
		if (    ( $filterString =~ /PASS/ )
			and ( $nFT =~ /LowQual/ )
			and ( $svlen >= $lengthOfSV ) )
		{
			$filterFlag = 2;
			$dupFlag    = 2;
		}
		if (
			    ( $svlen >= $lengthOfSV )
			and ( $mapq >= $overallMAPQ_Hotspot )
			and ( $peReads >= $overallSupportingReadsPE_Hotspot )
			and ( $srReads >= $overallSupportingReadsSR_Hotspot )

			and ( $tDV >= $sampleTumorSupportingReadsPE_Hotspot )
			and ( $nDV <= $sampleNormalSupportingReadsPE_Hotspot )
			#and ( $tGQ >= $sampleTumorGenotypeQualityFilterHotspot )
			#and ( $nGQ <= $sampleNormalGenotypeQualityFilter )
			and ( $dupFlag == 1 )
		  )
		{
			$filterFlag = 2;
		}
	}
	else {
		if (
			    ( $svlen >= $lengthOfSV )
			and ( $mapq >= $overallMAPQ )
			and ( $peReads >= $overallSupportingReadsPE )
			and ( $srReads >= $overallSupportingReadsSR )

			#and ( $svTNratio >= $tnRatio )
			and ( $tDV >= $sampleTumorSupportingReadsPE )
			and ( $nDV <= $sampleNormalSupportingReadsPE )
			#nd ( $tGQ >= $sampleTumorGenotypeQualityFilter )
			#and ( $nGQ <= $sampleNormalGenotypeQualityFilter )
		  )
		{
			$filterFlag = 2;
		}
	}
	return ($filterFlag);
}
#####################################
#####################################
#Filter record for translocations
sub FilterTranslocation {
	my (
		$svChr, $svStart, $svChr2,  $svEnd,   $filterString,
		$svlen, $mapq,    $peReads, $srReads, $tFT,
		$nFT,   $tGQ,     $nGQ,     $tDV,     $nDV
	) = @_;
	my $filterFlag = 1;

	#Start Checking if its Hotspot or not.
	my $hotspotFlag = CheckHotspot( $svChr, $svStart, $svChr2, $svEnd );
	if ( $hotspotFlag == 2 ) {

		#$svTNratio = $tnRatioHotspot if (! $svTNratio);
		my $dupFlag = 1;
		if (    ( $filterString =~ /PASS/ )
			and ( $nFT =~ /LowQual/ ) )
		{
			$filterFlag = 2;
			$dupFlag    = 2;
		}
		if (
			    ( $mapq >= $overallMAPQ_Hotspot )
			and ( $peReads >= $overallSupportingReadsPE_Hotspot )
			and ( $srReads >= $overallSupportingReadsSR_Hotspot )
			and ( $tDV >= $sampleTumorSupportingReadsPE_Hotspot )
			and ( $nDV <= $sampleNormalSupportingReadsPE_Hotspot )

			#and ( $svTNratio >= $tnRatioHotspot )
			#and ( $tGQ >= $sampleTumorGenotypeQualityFilterHotspot )
			#and ( $nGQ <= $sampleNormalGenotypeQualityFilter )
			and ( $dupFlag == 1 )
		  )
		{
			$filterFlag = 2;
		}
	}
	else {
		if (
			    ( $mapq >= $overallMAPQ )
			and ( $peReads >= $overallSupportingReadsPE )
			and ( $srReads >= $overallSupportingReadsSR )

			#and ( $svTNratio >= $tnRatio )
			and ( $tDV >= $sampleTumorSupportingReadsPE )
			and ( $nDV <= $sampleNormalSupportingReadsPE )
			#and ( $tGQ >= $sampleTumorGenotypeQualityFilter )
			#and ( $nGQ <= $sampleNormalGenotypeQualityFilter )
		  )
		{
			$filterFlag = 2;
		}
	}
	return ($filterFlag);
}
