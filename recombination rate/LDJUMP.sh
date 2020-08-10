#!/bin/bash -l
#SBATCH -A snic2020-5-20
#SBATCH -p core -n 5
#SBATCH -J LD_Jump
#SBATCH -t 30:00:00
#SBATCH -o LD_Jump.out
#SBATCH -e LD_Jump.err
#SBATCH --mail-user venkat.talla@ebc.uu.se
#SBATCH --mail-type=ALL


module load R/3.6.0 R_packages/3.6.0

cd /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts
Rscript LDJUMP.R

