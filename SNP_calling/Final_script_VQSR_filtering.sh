#!/bin/bash -l
#SBATCH -A snic2019-3-35
#SBATCH -p core -n 20
#SBATCH -J final
#SBATCH -t 100:00:00
#SBATCH -o Final_vcf.out
#SBATCH -e Final_vcf.err
#SBATCH --mail-user venkat.talla@ebc.uu.se
#SBATCH --mail-type=ALL


module load bioinfo-tools
module load samtools/1.3
module load bwa/0.7.12
module add FastQC/0.11.2
module add cutadapt/1.8.0
module add TrimGalore/0.4.0
module add vcftools

#java -Xmx120g -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/stm146.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/NJ116.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/mex919.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/NJ203.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/mex915.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/T14.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/mex536.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/mex986.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/NJ1.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/T9.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/PL10.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/PL7.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/PL6.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/PL4.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/PL5.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/PL1.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/PL8.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/PL9.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/PL2.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/HH9.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/HH5.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/HH6.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/HH7.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/HH8.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/HH10.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/HH3.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/HH2.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/HH4.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/HH1.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/ESB7.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/ESB6.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/ESB5.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/ESB9.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/ESB8.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/ESB1.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/ESB4.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/ESB2.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/ESB3.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/stm163.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/mex1527.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/HI033.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/HI023.final.bam.g.vcf \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/ESB10.final.bam.g.vcf \
#-o /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/genotype_all_20190315.vcf
#
#java -Xmx128g -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar \
#-R /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa \
#-T SelectVariants \
#-V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/genotype_all_20190315.vcf \
#-o /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/raw_indels_all_20190315.vcf \
#-selectType INDEL
#
#java -Xmx128g -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar \
#-R /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa \
#-T SelectVariants -V /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/genotype_all_20190315.vcf \
#-o /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/raw_snps_all_20190315.vcf \
#-selectType SNP
#
#minQ_indels=`perl /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/pl/find_qual.pl /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/raw_indels_all_20190315.vcf`
#minQ_snps=`perl /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/pl/find_qual.pl /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/raw_snps_all_20190315.vcf`
#
#vcftools --vcf /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/raw_indels_all_20190315.vcf --minQ $minQ_indels --recode --out /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/raw_indels_all_20190315_qual_filt
#vcftools --vcf /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/raw_snps_all_20190315.vcf --minQ $minQ_snps --recode --out /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/raw_snps_all_20190315_qual_filt                                           
#
java -Xmx128g -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa \
-nt 15 \
-mode SNP \
-input /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/raw_snps_all_20190315.vcf \
-resource:known=false,training=true,truth=true,prior=15.0 /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/raw_snps_all_20190315_qual_filt.recode.vcf \
-recalFile /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/recal/snp.tranches.recal \
-tranchesFile /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/recal/snp.tranches \
-an QD \
-an MQ \
-an MQRankSum \
-an ReadPosRankSum \
-an FS \
-an DP \
-tranche 100.0 \
-tranche 99.9 \
-tranche 99.0 \
-tranche 90.0 \
-rscriptFile /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/recal/snp.vqsr.r \
-allPoly

java -Xmx128g -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa \
-input /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/raw_snps_all_20190315.vcf \
-recalFile /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/recal/snp.tranches.recal \
-tranchesFile /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/recal/snp.tranches \
-o /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/recal_snps_all_20190315.vcf \
-mode SNP

java -Xmx128g -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar -T VariantRecalibrator \
-R /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa \
-nt 4 \
-mode INDEL \
-mG 4 \
-std 10.0 \
-input /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/raw_indels_all_20190315.vcf \
-resource:known=false,training=true,truth=true,prior=12.0 /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/raw_indels_all_20190315_qual_filt.recode.vcf  \
-recalFile /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/recal/indel.tranches.recal \
-tranchesFile /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/recal/indel.tranches \
-an QD \
-an DP \
-an FS \
-an ReadPosRankSum \
-an MQRankSum \
-tranche 100.0 \
-tranche 99.9 \
-tranche 99.9 \
-tranche 90.0 \
-rscriptFile /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/recal/indel.vqsr.r \
-allPoly

java -Xmx128g -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa \
-input /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/raw_indels_all_20190315.vcf \
-recalFile /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/recal/indel.tranches.recal \
-tranchesFile /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/recal/indel.tranches \
-o /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/recal_indels_all_20190315.vcf \
-mode INDEL

echo Job is done
date

