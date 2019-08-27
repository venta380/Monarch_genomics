#!/bin/bash -l
#SBATCH -A snic2019-3-35
#SBATCH -p node -n 1
#SBATCH -J Gvcfs-stm146
#SBATCH -t 60:00:00
#SBATCH -o Gvcfs-stm146.out
#SBATCH -e Gvcfs-stm146.err
#SBATCH --mail-user venkat.talla@ebc.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load samtools/1.2
module load bwa/0.7.12
module load bcftools
module load vcftools
module load java
module load cutadapt
ulimit -c unlimited
#java -Xmx32g -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar -T BaseRecalibrator -I /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/stm146.realign.bam -R /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa -o /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/stm146.bam.dedup.realign.calibration.csv -knownSites /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/golden_set.vcf.recode.vcf
#java  -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar -T PrintReads -I /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/stm146.realign.bam -R /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa -BQSR /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/stm146.bam.dedup.realign.calibration.csv -o /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/stm146.final.bam
#samtools flagstat /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/stm146.final.bam > /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/stm146.flagstat
java  -jar  /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar -T HaplotypeCaller -R /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa -nct 20 -I /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/stm146.final.bam --emitRefConfidence GVCF  -o /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/stm146.final.bam.g.vcf

