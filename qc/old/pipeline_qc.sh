#!/bin/sh

sm=$1
gp=$2
pr=$3
printf "\n"

printf "ALIGNMENT CHECKS:\n"
## check that solexa_PE.pl ran to completion for every sample/run combo
printf "Number of samples (including duplicates from different runs)                                              \t"
if ls $sm &> /dev/null; then cat $sm | wc -l; else printf "0\n"; fi
printf "  Number of filtered bams (*/*/*/*_FILTERED.bam)                                                           \t"
if ls */*/*/*_FILTERED.bam &> /dev/null; then ls */*/*/*_FILTERED.bam | wc -l; else printf "0\n"; fi
printf "  Number of filtered bam indexes (*/*/*/*FILTERED.bai)                                                    \t"
if ls */*/*/*_FILTERED.bai &> /dev/null; then ls */*/*/*_FILTERED.bai | wc -l; else printf "0\n"; fi

printf "\n"

## check that multiple samples from different runs were merged
printf "Number of unique samples                                                                                  \t"
if ls $sm &> /dev/null; then cut -f 2 $sm | sort | uniq | wc -l; else printf "0\n"; fi
printf "  Number of MD bams (*/*/*MD.bam)                                                                         \t"
if ls */*/*MD.bam &> /dev/null; then ls */*/*MD.bam | wc -l; else printf "0\n"; fi
printf "  Number of MD bam indexes (*/*/*MD.bai)                                                                  \t"
if ls */*/*MD.bai &> /dev/null; then ls */*/*MD.bai | wc -l; else printf "0\n"; fi

## check that bams were separated by sample after indel recalibration 
printf "  Number of indelRealigned recalibrated sample level bams (*indelRealigned_recal_*.bam)                   \t"
if ls *indelRealigned_recal_*.bam &> /dev/null; then ls *indelRealigned_recal_*.bam | wc -l; else printf "0\n"; fi
printf "  Number of indelRealigned recalibrated sample level bam indexes (*indelRealigned_recal_*.bai)            \t"
if ls *indelRealigned_recal_*.bai &> /dev/null; then ls *indelRealigned_recal_*.bai | wc -l; else printf "0\n"; fi

printf "\n"

## check that indel recalibration ran to completion
printf "Number of sample groups                                                                                   \t"
if ls $gp &> /dev/null; then cut -f 2 $gp | sort | uniq | wc -l; else printf "0\n"; fi
printf "  Number of grouped indel realigned recalibrated bams (*indelRealigned_recal.bam | grep -v CHR)           \t"
if ls *indelRealigned_recal.bam &> /dev/null; then ls *indelRealigned_recal.bam | grep -v CHR | wc -l; else printf "0\n"; fi
printf "  Number of grouped indel realigned recalibrated bam indexes (*indelRealigned_recal.bam | grep -v CHR)    \t"
if ls *indelRealigned_recal.bai &> /dev/null; then ls *indelRealigned_recal.bai | grep -v CHR | wc -l; else printf "0\n"; fi

printf "\n"

if [[ -z $gp ]] 
  then
    printf "No grouping file given. Not checking for variation detection results.\n"
  else
    printf "VARIATION DETECTION CHECKS:\n"
    ## check that UnifiedGenotyper and downstream recalibration/filtering ran to completion
    printf "UnifiedGenotyper files\n"
    printf "  *UnifiedGenotyper_RAW.vcf                                                                               \t"
    if ls *UnifiedGenotyper_RAW.vcf &> /dev/null; then printf "GOOD\n"; else printf "DOES NOT EXIST\n"; fi
    printf "  *UnifiedGenotyper_SNP_vqsr.vcf                                                                          \t"
    if ls *UnifiedGenotyper_SNP_vqsr.vcf &> /dev/null; then printf "GOOD\n"; else printf "DOES NOT EXIST\n"; fi
    printf "  *UnifiedGenotyper_vqsr.vcf                                                                              \t"
    if ls *UnifiedGenotyper_vqsr.vcf &> /dev/null; then printf "GOOD\n"; else printf "DOES NOT EXIST\n"; fi
    printf "  *UnifiedGenotyper.vcf                                                                                   \t"
    if ls *UnifiedGenotyper.vcf &> /dev/null; then printf "GOOD\n"; else printf "DOES NOT EXIST\n"; fi
    printf "  *UnifiedGenotyper.vcf_PAIRED_TCGA_MAF*                                                                  \t"
    if ls *UnifiedGenotyper.vcf_PAIRED_TCGA_MAF* &> /dev/null; then ls *UnifiedGenotyper.vcf_PAIRED_TCGA_MAF* | wc -l; else printf "0\n"; fi                                      
    printf "  *UnifiedGenotyper.vcf_UNPAIRED_TCGA_MAF*                                                                \t"                                     
    if ls *UnifiedGenotyper.vcf_UNPAIRED_TCGA_MAF* &> /dev/null; then ls *UnifiedGenotyper.vcf_UNPAIRED_TCGA_MAF* | wc -l; else printf "0\n"; fi   
    printf "\n" 

    ## check that HaplotypeCaller and downstream recalibration/filtering ran to completion
    printf "HaplotypeCaller files\n"
    printf "  *HaplotypeCaller_RAW.vcf                                                                                \t"
    if ls *HaplotypeCaller_RAW.vcf &> /dev/null; then printf "GOOD\n"; else printf "DOES NOT EXIST\n"; fi
    printf "  *HaplotypeCaller_SNP_vqsr.vcf                                                                           \t"
    if ls *HaplotypeCaller_SNP_vqsr.vcf &> /dev/null; then printf "GOOD\n"; else printf "DOES NOT EXIST\n"; fi
    printf "  *HaplotypeCaller_vqsr.vcf                                                                               \t"
    if ls *HaplotypeCaller_vqsr.vcf &> /dev/null; then printf "GOOD\n"; else printf "DOES NOT EXIST\n"; fi
    printf "  *HaplotypeCaller.vcf                                                                                    \t"
    if ls *HaplotypeCaller.vcf &> /dev/null; then printf "GOOD\n"; else printf "DOES NOT EXIST\n"; fi
    printf "  *HaplotypeCaller.vcf_PAIRED_TCGA_MAF*                                                                   \t"                                     
    if ls *HaplotypeCaller.vcf_PAIRED_TCGA_MAF* &> /dev/null; then ls *HaplotypeCaller.vcf_PAIRED_TCGA_MAF* | wc -l; else printf "0\n"; fi                 
    printf "  *HaplotypeCaller.vcf_UNPAIRED_TCGA_MAF*                                                                 \t"
    if ls *HaplotypeCaller.vcf_UNPAIRED_TCGA_MAF* &> /dev/null; then ls *HaplotypeCaller.vcf_UNPAIRED_TCGA_MAF* | wc -l; else printf "0\n"; fi  
fi

printf "\n"

if [[ -z "$pr" ]] 
  then
    printf "No pairing file given. Not checking for somatic analysis results.\n"
  else
    printf "SOMATIC ANALYSIS CHECKS:\n"
    printf "Number of sample pairs                                                                                    \t"
    if ls $pr &> /dev/null; then cat $pr | grep -v -i -e "NA$" -i -e "^NA" | wc -l; else printf "0\n"; fi
    printf "  Number of mutect files (should equal number of pairs)\n"
    printf "    *mutect_calls.txt                                                                                     \t"
    if ls *mutect_calls.txt &> /dev/null; then ls *mutect_calls.txt | wc -l; else printf "0\n"; fi
    printf "    *mutect_calls.vcf                                                                                     \t"
    if ls *mutect_calls.vcf &> /dev/null; then ls *mutect_calls.vcf | wc -l; else printf "0\n"; fi
    printf "    *mutect_calls.idx                                                                                     \t"
    if ls *mutect_calls.vcf.idx &> /dev/null; then ls *mutect_calls.vcf.idx | wc -l; else printf "0\n"; fi
    printf "  Number of varscan files (should equal 2*number of pairs)\n"
    printf "    *varscan*indel.vcf                                                                                    \t"
    if ls *varscan*indel* &> /dev/null; then ls *varscan*indel* | wc -l; else printf "0\n"; fi
    printf "    *varscan*snp.vcf                                                                                      \t"
    if ls *varscan*snp* &> /dev/null; then ls *varscan*snp* | wc -l; else printf "0\n"; fi
    printf "  Number of somatic sniper files (should equal 2*number of pairs)\n"
    printf "    *somatic_sniper.vcf                                                                                   \t"
    if ls *somatic_sniper.vcf &> /dev/null; then ls *somatic_sniper.vcf | wc -l; else printf "0\n"; fi
    printf "    *somatic_sniper_MAF.txt                                                                               \t"
    if ls *somatic_sniper_MAF.txt &> /dev/null; then ls *somatic_sniper_MAF.txt | wc -l; else printf "0\n"; fi
    printf "  Number of strelka files (should equal 4*number of pairs)\n"
    printf "    strelka/*/results/*.vcf                                                                               \t"
    if ls strelka/*/results/*.vcf &> /dev/null; then ls strelka/*/results/*.vcf | wc -l; else printf "0\n"; fi
    printf "    *STRELKA*MAF.txt                                                                                      \t"
    if ls *STRELKA*MAF.txt &> /dev/null; then ls *STRELKA*MAF.txt | wc -l; else printf "0\n"; fi
    printf "  Number of scalpel files (should equal number of pairs)\n"
    printf "    scalpel/*/somatic*indel.txt                                                                           \t"
    if ls scalpel/*/somatic*.indel.txt &> /dev/null; then ls scalpel/*/somatic*.indel.txt | wc -l; else printf "0\n"; fi
    #printf "    *scalpel_MAF.txt                                                                                      \t"
    #if ls *scalpel_MAF.txt &> /dev/null; then ls *scalpel_MAf.txt | wc -l; else printf "0\n"; fi
    printf "  Number of virmid files (should equal 2*number of pairs)\n"
    printf "    virmid/*/*.som.*.vcf                                                                                  \t"
    if ls virmid/*/*.som.*.vcf &> /dev/null; then ls virmid/*/*.som.*.vcf | wc -l; else printf "0\n"; fi    
    #printf "    *virmid_MAF.txt                                                                                       \t"
    #if ls *virmid_MAF.txt &> /dev/null; then ls *virmid_MAF.txt | wc -l; else printf "0\n"; fi
fi

printf "\n"

