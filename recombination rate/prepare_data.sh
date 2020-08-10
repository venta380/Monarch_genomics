#!/bin/bash -l
#SBATCH -A snic2020-5-20
#SBATCH -p core -n 20
#SBATCH -J LDhelmet
#SBATCH -t 200:00:00
#SBATCH -o LDhelmet.out
#SBATCH -e LDhelmet.err

module load bioinfo-tools
module load vcftools 
module load bcftools
module load tabix
module load LDhelmet

module load bioinfo-tools
module load vcftools 

cat /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa | vcf-consensus good_pos.recode.vcf -H 1  -s stm163 > stm163_1.fa
cat /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa | vcf-consensus good_pos.recode.vcf -H 1  -s stm146 > stm146_1.fa
cat /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa | vcf-consensus good_pos.recode.vcf -H 1  -s T9 > T9_1.fa
cat /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa | vcf-consensus good_pos.recode.vcf -H 1  -s T14 > T14_1.fa
cat /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa | vcf-consensus good_pos.recode.vcf -H 1  -s NJ203 > NJ203_1.fa
cat /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa | vcf-consensus good_pos.recode.vcf -H 1  -s NJ116 > NJ116_1.fa
cat /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa | vcf-consensus good_pos.recode.vcf -H 1  -s NJ1 > NJ1_1.fa
cat /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa | vcf-consensus good_pos.recode.vcf -H 1  -s HI023 > HI023_1.fa
cat /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa | vcf-consensus good_pos.recode.vcf -H 1  -s HI033 > HI033_1.fa
cat /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa | vcf-consensus good_pos.recode.vcf -H 1  -s mex986 > mex986_1.fa
cat /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa | vcf-consensus good_pos.recode.vcf -H 1  -s mex919 > mex919_1.fa
cat /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa | vcf-consensus good_pos.recode.vcf -H 1  -s mex915 > mex915_1.fa
cat /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa | vcf-consensus good_pos.recode.vcf -H 1  -s mex536 > mex536_1.fa
cat /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa | vcf-consensus good_pos.recode.vcf -H 1  -s mex1527 > mex1527_1.fa

cat /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa | vcf-consensus good_pos.recode.vcf -H 2  -s stm163 > stm163_2.fa
cat /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa | vcf-consensus good_pos.recode.vcf -H 2  -s stm146 > stm146_2.fa
cat /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa | vcf-consensus good_pos.recode.vcf -H 2  -s T9 > T9_2.fa
cat /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa | vcf-consensus good_pos.recode.vcf -H 2  -s T14 > T14_2.fa
cat /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa | vcf-consensus good_pos.recode.vcf -H 2  -s NJ203 > NJ203_2.fa
cat /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa | vcf-consensus good_pos.recode.vcf -H 2  -s NJ116 > NJ116_2.fa
cat /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa | vcf-consensus good_pos.recode.vcf -H 2  -s NJ1 > NJ1_2.fa
cat /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa | vcf-consensus good_pos.recode.vcf -H 2  -s HI023 > HI023_2.fa
cat /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa | vcf-consensus good_pos.recode.vcf -H 2  -s HI033 > HI033_2.fa
cat /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa | vcf-consensus good_pos.recode.vcf -H 2  -s mex986 > mex986_2.fa
cat /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa | vcf-consensus good_pos.recode.vcf -H 2  -s mex919 > mex919_2.fa
cat /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa | vcf-consensus good_pos.recode.vcf -H 2  -s mex915 > mex915_2.fa
cat /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa | vcf-consensus good_pos.recode.vcf -H 2  -s mex536 > mex536_2.fa
cat /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa | vcf-consensus good_pos.recode.vcf -H 2  -s mex1527 > mex1527_2.fa

python Generate_data_for_LD_helmet.py

