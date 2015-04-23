
# Get Targeted Mutations, generate OLD style MAF with pairing info
# Filter for somatic events and remove LowQual events
/home/mpirun/bin/bedtools intersect -a Proj_3213/Proj_3213_COMBINED_VARIANTS.vcf -b /ifs/data/mpirun/data/targets/hg19__MegaGene__v120309.bed -header | /home/mpirun/tools/maf/vcf2mafV3.py -p pairing.txt | /home/mpirun/tools/maf/pA_qSom0+noLow.py | /usr/bin/tee maf0.txt

# Convert old (NDS) MAF to official TCGA MAF
## NOTE; DON'T FORGET <> AROUND INPUT FILE (maf0.txt)
/home/mpirun/tools/maf/oldMAF2tcgaMAF.py <maf0.txt> maf1.txt

# Annotate with Oncotator
## expects db.properties to be in working directory
ln -s /home/mpirun/bin/oncotator/db.properties .
/home/mpirun/bin/oncotator/oncotateMaf.sh maf1.txt maf2.txt

# Create annotated TCGA MAF
cat maf2.txt |/home/mpirun/tools/maf/pA_Functional.py | /home/mpirun/tools/maf/pA_fixHugo.py >jabobl__p3213__NoLow,qSom0,Functional___TCGA_MAF.txt

# Create nice MAF with essential columns
cat maf2.txt | /home/mpirun/tools/maf/pA_Functional.py | /home/mpirun/tools/maf/pA_reSortCols.txt /home/mpirun/tools/finalCols.txt >jabobl__p3213__NoLow,qSom0,Functional___MAF3.txt

