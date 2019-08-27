import sys
import os
import string
import pandas
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import seaborn as sns
import itertools
import math
import time
import sys
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
import gzip
from subprocess import Popen, PIPE




def run_BWA_GATK_lib_2(Sample_ID, genome, out_dir, fastq_1, fastq_2, fastq_3, fastq_4, outfile):
	outFileName = outfile
	outFile = open(outFileName, "w")
	outFile.write('#!/bin/bash -l'+'\n')
	outFile.write('#SBATCH -A snic2019-3-35'+'\n')
	outFile.write('#SBATCH -p core -n 6'+'\n')
	outFile.write('#SBATCH -J mapping-'+Sample_ID+'\n')
	outFile.write('#SBATCH -t 30:00:00'+'\n')
	outFile.write('#SBATCH -o mapping-'+Sample_ID+'.out'+'\n')
	outFile.write('#SBATCH -e mapping-'+Sample_ID+'.err'+'\n')
	outFile.write('#SBATCH --mail-user venkat.talla@ebc.uu.se'+'\n')
	outFile.write('#SBATCH --mail-type=ALL'+'\n')
	outFile.write('\n')
	outFile.write("module load bioinfo-tools"+'\n')
	outFile.write("module load samtools/1.2"+'\n')
	outFile.write("module load bwa/0.7.12"+'\n')
	outFile.write("module load bcftools"+'\n')
	outFile.write("module load vcftools"+'\n')
	outFile.write("module load java"+'\n')
	outFile.write("module load cutadapt"+'\n')
	outFile.write("ulimit -c unlimited"+'\n')
	outFile.write('cd $TMPDIR')
	outFile.write('\n')
	outFile.write('wget '+fastq_1[0]+' -O '+fastq_1[1]+' &'+'\n')
	outFile.write('wget '+fastq_2[0]+' -O '+fastq_2[1]+' &'+'\n')
	outFile.write('wget '+fastq_1[0]+' -O '+fastq_3[1]+' &'+'\n')
	outFile.write('wget '+fastq_2[0]+' -O '+fastq_4[1]+' &'+'\n')
	outFile.write('wait'+'\n')
	outFile.write('\n')
	outFile.write("cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o "+'.'.join(fastq_1[1].split('.')[0:2])+'_trimmed.fastq.gz -p '+'.'.join(fastq_2[1].split('.')[0:2])+'_trimmed.fastq.gz '+' '+fastq_1[1]+' '+fastq_2[1]+' '+'\n')
	outFile.write("cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o "+'.'.join(fastq_3[1].split('.')[0:2])+'_trimmed.fastq.gz -p '+'.'.join(fastq_4[1].split('.')[0:2])+'_trimmed.fastq.gz '+' '+fastq_3[1]+' '+fastq_4[1]+' '+'\n')
	outFile.write('\n')
	outFile.write('\n')
	outFile.write('bwa mem -t 6 -M '+genome+' -R '+' '+'"@RG\\tLB:Lib1\\tID:1\\tSM:'+Sample_ID+'\\tPL:ILLUMINA"'+' $TMPDIR/'+'.'.join(fastq_1[1].split('.')[0:2])+'_trimmed.fastq.gz '+' $TMPDIR/'+'.'.join(fastq_2[1].split('.')[0:2])+'_trimmed.fastq.gz '+' '+'| samtools import /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa.fai - - | samtools sort - '+'$TMPDIR/'+Sample_ID+'_1'+'\n')
	outFile.write('samtools index '+'$TMPDIR/'+Sample_ID+'_1'+'.bam'+'\n')
	outFile.write('bwa mem -t 6 -M '+genome+' -R '+' '+'"@RG\\tLB:Lib1\\tID:1\\tSM:'+Sample_ID+'\\tPL:ILLUMINA"'+' $TMPDIR/'+'.'.join(fastq_3[1].split('.')[0:2])+'_trimmed.fastq.gz '+' $TMPDIR/'+'.'.join(fastq_4[1].split('.')[0:2])+'_trimmed.fastq.gz '+' '+'| samtools import /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa.fai - - | samtools sort - '+'$TMPDIR/'+Sample_ID+'_2'+'\n')
	outFile.write('samtools index '+'$TMPDIR/'+Sample_ID+'_2'+'.bam'+'\n')
	outFile.write('java -Xmx32g -jar /sw/apps/bioinfo/picard/1.127/milou/picard.jar MarkDuplicates INPUT='+'$TMPDIR/'+Sample_ID+'_1'+'.bam '+'OUTPUT='+'$TMPDIR/'+Sample_ID+'_1'+'.bam.dedup.bam METRICS_FILE='+'$TMPDIR/'+Sample_ID+'_1'+'.bam.metrics REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT'+'\n')
	outFile.write('java -Xmx32g -jar /sw/apps/bioinfo/picard/1.127/milou/picard.jar MarkDuplicates INPUT='+'$TMPDIR/'+Sample_ID+'_2'+'.bam '+'OUTPUT='+'$TMPDIR/'+Sample_ID+'_2'+'.bam.dedup.bam METRICS_FILE='+'$TMPDIR/'+Sample_ID+'_2'+'.bam.metrics REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT'+'\n')
	outFile.write('samtools merge -cf '+'$TMPDIR/'+Sample_ID+'.bam.dedup.bam '+'$TMPDIR/'+Sample_ID+'_1'+'.bam.dedup.bam '+'$TMPDIR/'+Sample_ID+'_2'+'.bam.dedup.bam'+'\n')
	outFile.write('samtools index '+'$TMPDIR/'+Sample_ID+'.bam.dedup.bam'+'\n')
	outFile.write('java -Xmx32g -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar -I ' + '$TMPDIR/'+Sample_ID+'.bam.dedup.bam -R /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa -T RealignerTargetCreator -o '+'$TMPDIR/'+Sample_ID+'.intervals'+'\n')
	outFile.write('java -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar -I ' + '$TMPDIR/'+Sample_ID+'.bam.dedup.bam -R /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa -T IndelRealigner --filter_bases_not_stored -o '+out_dir+Sample_ID+'.realign.bam -targetIntervals '+'$TMPDIR/'+Sample_ID+'.intervals'+'\n')
	outFile.write('samtools index '+out_dir+Sample_ID+'.realign.bam'+'\n')
	outFile.write('\n')
	outFile.write('\n')


def run_BWA_GATK_lib_1(Sample_ID, genome, out_dir, fastq_1, fastq_2, outfile):
	outFileName = outfile
	outFile = open(outFileName, "w")
	outFile.write('#!/bin/bash -l'+'\n')
	outFile.write('#SBATCH -A snic2019-3-35'+'\n')
	outFile.write('#SBATCH -p core -n 10'+'\n')
	outFile.write('#SBATCH -J mapping-'+Sample_ID+'\n')
	outFile.write('#SBATCH -t 30:00:00'+'\n')
	outFile.write('#SBATCH -o mapping-'+Sample_ID+'.out'+'\n')
	outFile.write('#SBATCH -e mapping-'+Sample_ID+'.err'+'\n')
	outFile.write('#SBATCH --mail-user venkat.talla@ebc.uu.se'+'\n')
	outFile.write('#SBATCH --mail-type=ALL'+'\n')
	outFile.write('\n')
	outFile.write("module load bioinfo-tools"+'\n')
	outFile.write("module load samtools/1.2"+'\n')
	outFile.write("module load bwa/0.7.12"+'\n')
	outFile.write("module load bcftools"+'\n')
	outFile.write("module load vcftools"+'\n')
	outFile.write("module load java"+'\n')
	outFile.write("module load cutadapt"+'\n')
	outFile.write("ulimit -c unlimited"+'\n')
	outFile.write('cd $TMPDIR')
	outFile.write('\n')
	outFile.write('wget '+fastq_1[0]+' -O '+fastq_1[1]+' &'+'\n')
	outFile.write('wget '+fastq_2[0]+' -O '+fastq_2[1]+' &'+'\n')
	outFile.write('wait'+'\n')
	outFile.write('\n')
	outFile.write("cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o "+'.'.join(fastq_1[1].split('.')[0:2])+'_trimmed.fastq.gz -p '+'.'.join(fastq_2[1].split('.')[0:2])+'_trimmed.fastq.gz '+' '+fastq_1[1]+' '+fastq_2[1]+' '+'\n')
	outFile.write('\n')
	outFile.write('\n')
	outFile.write('bwa mem -t 10 -M '+genome+' -R '+' '+'"@RG\\tLB:Lib1\\tID:1\\tSM:'+Sample_ID+'\\tPL:ILLUMINA"'+' $TMPDIR/'+'.'.join(fastq_1[1].split('.')[0:2])+'_trimmed.fastq.gz '+' $TMPDIR/'+'.'.join(fastq_2[1].split('.')[0:2])+'_trimmed.fastq.gz '+' '+'| samtools import /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa.fai - - | samtools sort - '+'$TMPDIR/'+Sample_ID+'\n')
	outFile.write('samtools index '+'$TMPDIR/'+Sample_ID+'.bam'+'\n')
	outFile.write('java -Xmx32g -jar /sw/apps/bioinfo/picard/1.127/milou/picard.jar MarkDuplicates INPUT='+'$TMPDIR/'+Sample_ID+'.bam '+'OUTPUT='+'$TMPDIR/'+Sample_ID+'.bam.dedup.bam METRICS_FILE='+'$TMPDIR/'+Sample_ID+'.bam.metrics REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT'+'\n')
	outFile.write('samtools index '+'$TMPDIR/'+Sample_ID+'.bam.dedup.bam'+'\n')
	outFile.write('java -Xmx32g -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar -I ' + '$TMPDIR/'+Sample_ID+'.bam.dedup.bam -R /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa -T RealignerTargetCreator -o '+'$TMPDIR/'+Sample_ID+'.intervals'+'\n')
	outFile.write('java -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar -I ' + '$TMPDIR/'+Sample_ID+'.bam.dedup.bam -R /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa -T IndelRealigner --filter_bases_not_stored -o '+out_dir+Sample_ID+'.realign.bam -targetIntervals '+'$TMPDIR/'+Sample_ID+'.intervals'+'\n')
	outFile.write('samtools index '+out_dir+Sample_ID+'.realign.bam'+'\n')
	outFile.write('\n')
	outFile.write('\n')



def run_BCF_tools_snp_calling(Sample_ID, genome, out_dir, outfile):
	outFileName = outfile
	outFile = open(outFileName, "w")
	outFile.write('#!/bin/bash -l'+'\n')
	outFile.write('#SBATCH -A snic2019-3-35'+'\n')
	outFile.write('#SBATCH -p core -n 3'+'\n')
	outFile.write('#SBATCH -J SNP_known_sites-'+Sample_ID+'\n')
	outFile.write('#SBATCH -t 20:00:00'+'\n')
	outFile.write('#SBATCH -o SNP_known_sites-'+Sample_ID+'.out'+'\n')
	outFile.write('#SBATCH -e SNP_known_sites-'+Sample_ID+'.err'+'\n')
	outFile.write('#SBATCH --mail-user venkat.talla@ebc.uu.se'+'\n')
	outFile.write('#SBATCH --mail-type=ALL'+'\n')
	outFile.write('\n')
	outFile.write("module load bioinfo-tools"+'\n')
	outFile.write("module load samtools/1.2"+'\n')
	outFile.write("module load bwa/0.7.12"+'\n')
	outFile.write("module load bcftools"+'\n')
	outFile.write("module load vcftools"+'\n')
	outFile.write("bcftools mpileup -f "+genome+' '+out_dir+Sample_ID+'.realign.bam | bcftools call -mv -Ob -o '+'$TMPDIR/'+Sample_ID+'_calls.bcf'+'\n')
	outFile.write('bcftools view -i '"'%QUAL>=80'"' '+'$TMPDIR/'+Sample_ID+'_calls.bcf > '+out_dir+Sample_ID+'_calls.vcf'+'\n')
	outFile.write('vcftools --vcf '+out_dir+Sample_ID+'_calls.vcf --remove-indels --exclude-bed /proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa.out.bed --recode --recode-INFO-all --out '+out_dir+Sample_ID+'_SNPs_only'+'\n')




def make_final_bam_and_call_gVCF(Sample_ID, genome, out_dir, outfile, Knownsites):
	outFileName = outfile
	outFile = open(outFileName, "w")
	outFile.write('#!/bin/bash -l'+'\n')
	outFile.write('#SBATCH -A snic2019-3-35'+'\n')
	outFile.write('#SBATCH -p core -n 1'+'\n')
	outFile.write('#SBATCH -J Gvcfs-'+Sample_ID+'\n')
	outFile.write('#SBATCH -t 60:00:00'+'\n')
	outFile.write('#SBATCH -o Gvcfs-'+Sample_ID+'.out'+'\n')
	outFile.write('#SBATCH -e Gvcfs-'+Sample_ID+'.err'+'\n')
	outFile.write('#SBATCH --mail-user venkat.talla@ebc.uu.se'+'\n')
	outFile.write('#SBATCH --mail-type=ALL'+'\n')
	outFile.write('\n')
	outFile.write("module load bioinfo-tools"+'\n')
	outFile.write("module load samtools/1.2"+'\n')
	outFile.write("module load bwa/0.7.12"+'\n')
	outFile.write("module load bcftools"+'\n')
	outFile.write("module load vcftools"+'\n')
	outFile.write("module load java"+'\n')
	outFile.write("module load cutadapt"+'\n')
	outFile.write("ulimit -c unlimited"+'\n')
	outFile.write("#java -Xmx32g -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar -T BaseRecalibrator -I "+out_dir+Sample_ID+'.realign.bam -R '+genome+' -o '+out_dir+Sample_ID+'.bam.dedup.realign.calibration.csv -knownSites '+out_dir+Knownsites+'.recode.vcf'+'\n')
	outFile.write("#java  -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar -T PrintReads -I "+out_dir+Sample_ID+'.realign.bam -R '+genome+' -BQSR '+out_dir+Sample_ID+'.bam.dedup.realign.calibration.csv -o '+out_dir+Sample_ID+'.final.bam'+'\n')
	outFile.write("#samtools flagstat "+out_dir+Sample_ID+'.final.bam > '+out_dir+Sample_ID+'.flagstat'+'\n')
	outFile.write("java  -jar  /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar -T HaplotypeCaller -R "+genome+' -I '+out_dir+Sample_ID+'.final.bam --emitRefConfidence GVCF  -o '+out_dir+Sample_ID+'.final.bam.g.vcf'+'\n'+'\n')


def run_final_bam_and_call_gVCF(samples, Knownsites):
	samples=samples
	Knownsites=Knownsites
	for i in samples:
		Sample_ID=str(i)
		out_dir='/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/'
		genome='/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa'
		outfile='/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/'+Sample_ID+'.sh'
		make_final_bam_and_call_gVCF(Sample_ID, genome, out_dir, outfile, Knownsites)



samples_2_libraries=set(['HH2','HH3','HH5','HH9','PL1','PL4','PL5','PL6','PL7','PL9','PL10','ESB1','ESB2','ESB5','ESB6','ESB7'])
samples_1_libraries=set(['HH1','HH4','HH6','HH7','HH8','HH10','PL2','PL3','PL8','ESB3','ESB4','ESB8','ESB9','ESB10','stm163','stm146','T9','T14','NJ203','NJ116','NJ1','HI023','HI033','mex986','mex919','mex915','mex536','mex1527'])

link_dict={'HH1_R1':['https://www.dropbox.com/s/xh0zjvlgikfr713/HH1_trimmed_R1.fastq.gz?dl=0','HH1_trimmed_R1.fastq.gz'], 'HH1_R2':['https://www.dropbox.com/s/7f0bu9aeud2a4pv/HH1_trimmed_R2.fastq.gz?dl=0','HH1_trimmed_R2.fastq.gz'], 'HH2_R1':['https://www.dropbox.com/s/nws7aml7btdm3qm/HH2_trimmed_R1.fastq.gz?dl=0','HH2_trimmed_R1.fastq.gz'], 'HH2_R2':['https://www.dropbox.com/s/y9wjmxibbdug3nf/HH2_trimmed_R2.fastq.gz?dl=0','HH2_trimmed_R2.fastq.gz'], 'HH3_R1':['https://www.dropbox.com/s/1wo3rghvil63wdh/HH3_trimmed_R1.fastq.gz?dl=0','HH3_trimmed_R1.fastq.gz'], 'HH3_R2':['https://www.dropbox.com/s/2ewrq3r69rxgk5j/HH3_trimmed_R2.fastq.gz?dl=0','HH3_trimmed_R2.fastq.gz'], 'HH4_R1':['https://www.dropbox.com/s/envwqk7oyf7003m/HH4_trimmed_R1.fastq.gz?dl=0','HH4_trimmed_R1.fastq.gz'], 'HH4_R2':['https://www.dropbox.com/s/2wzgcv73off1wvy/HH4_trimmed_R2.fastq.gz?dl=0','HH4_trimmed_R2.fastq.gz'], 'HH5_R1':['https://www.dropbox.com/s/j9y3j1p9fltgvv2/HH5_trimmed_R1.fastq.gz?dl=0','HH5_trimmed_R1.fastq.gz'], 'HH5_R2':['https://www.dropbox.com/s/g9lvoxho3wqeit1/HH5_trimmed_R2.fastq.gz?dl=0','HH5_trimmed_R2.fastq.gz'], 'HH6_R1':['https://www.dropbox.com/s/lkq0t79xb9mvo57/HH6_trimmed_R1.fastq.gz?dl=0','HH6_trimmed_R1.fastq.gz'], 'HH6_R2':['https://www.dropbox.com/s/up5nfq3ksmv1iaf/HH6_trimmed_R2.fastq.gz?dl=0','HH6_trimmed_R2.fastq.gz'], 'HH7_R1':['https://www.dropbox.com/s/dc1wcxvtk5ws8hg/HH7_trimmed_R1.fastq.gz?dl=0','HH7_trimmed_R1.fastq.gz'], 'HH7_R2':['https://www.dropbox.com/s/wdrc1kvzmbgq9ux/HH7_trimmed_R2.fastq.gz?dl=0','HH7_trimmed_R2.fastq.gz'], 'HH8_R1':['https://www.dropbox.com/s/m07xhduuiieke1e/HH8_trimmed_R1.fastq.gz?dl=0','HH8_trimmed_R1.fastq.gz'], 'HH8_R2':['https://www.dropbox.com/s/1wf8p9bxqxcrx9p/HH8_trimmed_R2.fastq.gz?dl=0','HH8_trimmed_R2.fastq.gz'], 'HH9_R1':['https://www.dropbox.com/s/a4d05s7cc9g4r0d/HH9_trimmed_R1.fastq.gz?dl=0','HH9_trimmed_R1.fastq.gz'], 'HH9_R2':['https://www.dropbox.com/s/1wn1f2l3eei4a87/HH9_trimmed_R2.fastq.gz?dl=0','HH9_trimmed_R2.fastq.gz'], 'HH10_R1':['https://www.dropbox.com/s/y73g2uobcizs0vo/HH10_trimmed_R1.fastq.gz?dl=0','HH10_trimmed_R1.fastq.gz'], 'HH10_R2':['https://www.dropbox.com/s/wm4h46u08wycj24/HH10_trimmed_R2.fastq.gz?dl=0','HH10_trimmed_R2.fastq.gz'], 'PL1_R1':['https://www.dropbox.com/s/uhfnxr5px4qe4ac/PL1_trimmed_R1.fastq.gz?dl=0','PL1_trimmed_R1.fastq.gz'], 'PL1_R2':['https://www.dropbox.com/s/ndqovf1ew4x5au0/PL1_trimmed_R2.fastq.gz?dl=0','PL1_trimmed_R2.fastq.gz'], 'PL2_R1':['https://www.dropbox.com/s/qtjuwh7f9cjplz2/PL2_trimmed_R1.fastq.gz?dl=0','PL2_trimmed_R1.fastq.gz'], 'PL2_R2':['https://www.dropbox.com/s/qgnw4dbsr088fnk/PL2_trimmed_R2.fastq.gz?dl=0','PL2_trimmed_R2.fastq.gz'], 'PL3_R1':['https://www.dropbox.com/s/8valtzb6zpvok2u/PL3_trimmed_R1.fastq.gz?dl=0','PL3_trimmed_R1.fastq.gz'], 'PL3_R2':['https://www.dropbox.com/s/8bbib4bbfhcvxxy/PL3_trimmed_R2.fastq.gz?dl=0','PL3_trimmed_R2.fastq.gz'], 'PL4_R1':['https://www.dropbox.com/s/nbbxdapkyafp8at/PL4_trimmed_R1.fastq.gz?dl=0','PL4_trimmed_R1.fastq.gz'], 'PL4_R2':['https://www.dropbox.com/s/rc2o62imv8mj4k1/PL4_trimmed_R2.fastq.gz?dl=0','PL4_trimmed_R2.fastq.gz'], 'PL5_R1':['https://www.dropbox.com/s/8a3wrvb3tbuinq0/PL5_trimmed_R1.fastq.gz?dl=0','PL5_trimmed_R1.fastq.gz'], 'PL5_R2':['https://www.dropbox.com/s/79ybggy9qe139er/PL5_trimmed_R2.fastq.gz?dl=0','PL5_trimmed_R2.fastq.gz'], 'PL6_R1':['https://www.dropbox.com/s/8mmg3n1kkguna8j/PL6_trimmed_R1.fastq.gz?dl=0','PL6_trimmed_R1.fastq.gz'], 'PL6_R2':['https://www.dropbox.com/s/ee3tak4dwgct6ok/PL6_trimmed_R2.fastq.gz?dl=0','PL6_trimmed_R2.fastq.gz'], 'PL7_R1':['https://www.dropbox.com/s/xl4utzifbbf57wg/PL7_trimmed_R1.fastq.gz?dl=0','PL7_trimmed_R1.fastq.gz'], 'PL7_R2':['https://www.dropbox.com/s/fn6brju2pz4ef4n/PL7_trimmed_R2.fastq.gz?dl=0','PL7_trimmed_R2.fastq.gz'], 'PL8_R1':['https://www.dropbox.com/s/cg3s31qjkmkg2hg/PL8_trimmed_R1.fastq.gz?dl=0','PL8_trimmed_R1.fastq.gz'], 'PL8_R2':['https://www.dropbox.com/s/p1zlp52wyrrdm09/PL8_trimmed_R2.fastq.gz?dl=0','PL8_trimmed_R2.fastq.gz'], 'PL9_R1':['https://www.dropbox.com/s/z3f2hulxpw5i1ek/PL9_trimmed_R1.fastq.gz?dl=0','PL9_trimmed_R1.fastq.gz'], 'PL9_R2':['https://www.dropbox.com/s/zia2fiylemidx2c/PL9_trimmed_R2.fastq.gz?dl=0','PL9_trimmed_R2.fastq.gz'], 'PL10_R1':['https://www.dropbox.com/s/ju1odcufnwz2s6h/PL10_trimmed_R1.fastq.gz?dl=0','PL10_trimmed_R1.fastq.gz'], 'PL10_R2':['https://www.dropbox.com/s/dsd5r7irbf79wcg/PL10_trimmed_R2.fastq.gz?dl=0','PL10_trimmed_R2.fastq.gz'], 'ESB1_R1':['https://www.dropbox.com/s/ut9gpj6t40a6u7f/ESB1_trimmed_R1.fastq.gz?dl=0','ESB1_trimmed_R1.fastq.gz'], 'ESB1_R2':['https://www.dropbox.com/s/btnt42rlmvarxxw/ESB1_trimmed_R2.fastq.gz?dl=0','ESB1_trimmed_R2.fastq.gz'], 'ESB2_R1':['https://www.dropbox.com/s/7v3mxedrz9j6jne/ESB2_trimmed_R1.fastq.gz?dl=0','ESB2_trimmed_R1.fastq.gz'], 'ESB2_R2':['https://www.dropbox.com/s/depq1ubjt201m1i/ESB2_trimmed_R2.fastq.gz?dl=0','ESB2_trimmed_R2.fastq.gz'], 'ESB3_R1':['https://www.dropbox.com/s/mhjdbxcsv5iepfw/ESB3_trimmed_R1.fastq.gz?dl=0','ESB3_trimmed_R1.fastq.gz'], 'ESB3_R2':['https://www.dropbox.com/s/qf8h85cciwoj4ze/ESB3_trimmed_R2.fastq.gz?dl=0','ESB3_trimmed_R2.fastq.gz'], 'ESB4_R1':['https://www.dropbox.com/s/xu9nmrgrzcmjxp1/ESB4_trimmed_R1.fastq.gz?dl=0','ESB4_trimmed_R1.fastq.gz'], 'ESB4_R2':['https://www.dropbox.com/s/mc1w7mt7qirk1a1/ESB4_trimmed_R2.fastq.gz?dl=0','ESB4_trimmed_R2.fastq.gz'], 'ESB5_R1':['https://www.dropbox.com/s/d3lkn3h0wvzrqum/ESB5_trimmed_R1.fastq.gz?dl=0','ESB5_trimmed_R1.fastq.gz'], 'ESB5_R2':['https://www.dropbox.com/s/u95nef5rf3kr636/ESB5_trimmed_R2.fastq.gz?dl=0','ESB5_trimmed_R2.fastq.gz'], 'ESB6_R1':['https://www.dropbox.com/s/whe8zq33utwpyz7/ESB6_trimmed_R1.fastq.gz?dl=0','ESB6_trimmed_R1.fastq.gz'], 'ESB6_R2':['https://www.dropbox.com/s/wwd1454w0lqwbtz/ESB6_trimmed_R2.fastq.gz?dl=0','ESB6_trimmed_R2.fastq.gz'], 'ESB7_R1':['https://www.dropbox.com/s/6cg9votvdbvfrlq/ESB7_trimmed_R1.fastq.gz?dl=0','ESB7_trimmed_R1.fastq.gz'], 'ESB7_R2':['https://www.dropbox.com/s/53jgup20y9e8t8o/ESB7_trimmed_R2.fastq.gz?dl=0','ESB7_trimmed_R2.fastq.gz'], 'ESB8_R1':['https://www.dropbox.com/s/gl90fofj2p149he/ESB8_trimmed_R1.fastq.gz?dl=0','ESB8_trimmed_R1.fastq.gz'], 'ESB8_R2':['https://www.dropbox.com/s/b7xlcmcfjuf1tnf/ESB8_trimmed_R2.fastq.gz?dl=0','ESB8_trimmed_R2.fastq.gz'], 'ESB9_R1':['https://www.dropbox.com/s/7t1q32306lmav35/ESB9_trimmed_R1.fastq.gz?dl=0','ESB9_trimmed_R1.fastq.gz'], 'ESB9_R2':['https://www.dropbox.com/s/dmfrakc19q4og6s/ESB9_trimmed_R2.fastq.gz?dl=0','ESB9_trimmed_R2.fastq.gz'], 'ESB10_R1':['https://www.dropbox.com/s/l4rnuv14bkvmp4v/ESB10_trimmed_R1.fastq.gz?dl=0','ESB10_trimmed_R1.fastq.gz'], 'ESB10_R2':['https://www.dropbox.com/s/ilrb62cjef46uv9/ESB10_trimmed_R2.fastq.gz?dl=0','ESB10_trimmed_R2.fastq.gz'], 'stm163_R1':['https://www.dropbox.com/s/8sfjy2wnvpjhija/stm163_trimmed_1.fastq.gz?dl=0', 'stm163_trimmed_1.fastq.gz'], 'stm163_R2':['https://www.dropbox.com/s/7kue6jtj2ize63e/stm163_trimmed_2.fastq.gz?dl=0', 'stm163_trimmed_2.fastq.gz'], 'stm146_R1':['https://www.dropbox.com/s/fhmj1furxx0qx33/stm146_trimmed_R1.fastq.gz?dl=0', 'stm146_trimmed_R1.fastq.gz'], 'stm146_R2':['https://www.dropbox.com/s/t6bvjxm7d3qe42g/stm146_trimmed_R2.fastq.gz?dl=0', 'stm146_trimmed_R2.fastq.gz'], 'T9_R1':['https://www.dropbox.com/s/lvc9fuwi5ftfvik/T9_trimmed_R1.fastq.gz?dl=0', 'T9_trimmed_R1.fastq.gz'], 'T9_R2':['https://www.dropbox.com/s/q487gymddf14k3z/T9_trimmed_R2.fastq.gz?dl=0', 'T9_trimmed_R2.fastq.gz'], 'T14_R1':['https://www.dropbox.com/s/5tn9iry58mxtgdh/T14_trimmed_R1.fastq.gz?dl=0', 'T14_trimmed_R1.fastq.gz'], 'T14_R2':['https://www.dropbox.com/s/woofqmcdsfrtopv/T14_trimmed_R2.fastq.gz?dl=0', 'T14_trimmed_R2.fastq.gz'], 'NJ203_R1':['https://www.dropbox.com/s/utx7dl9w3d3auq4/NJ203_trimmed_R1.fastq.gz?dl=0', 'NJ203_trimmed_R1.fastq.gz'], 'NJ203_R2':['https://www.dropbox.com/s/icq28x1kaxq95sp/NJ203_trimmed_R2.fastq.gz?dl=0', 'NJ203_trimmed_R2.fastq.gz'], 'NJ116_R1':['https://www.dropbox.com/s/djmbckc54b3p3kt/NJ116_trimmed_R1.fastq.gz?dl=0', 'NJ116_trimmed_R1.fastq.gz'], 'NJ116_R2':['https://www.dropbox.com/s/v6765pttzs7yujl/NJ116_trimmed_R2.fastq.gz?dl=0', 'NJ116_trimmed_R2.fastq.gz'], 'NJ1_R1':['https://www.dropbox.com/s/eq1qwgy4c69nx8f/NJ1_trimmed_R1.fastq.gz?dl=0', 'NJ1_trimmed_R1.fastq.gz'], 'NJ1_R2':['https://www.dropbox.com/s/88yo5nwqjoknmhr/NJ1_trimmed_R2.fastq.gz?dl=0', 'NJ1_trimmed_R2.fastq.gz'], 'HI023_R1':['https://www.dropbox.com/s/fg2w45bygq6kt6f/HI023_trimmed_1.fastq.gz?dl=0', 'HI023_trimmed_1.fastq.gz'], 'HI023_R2':['https://www.dropbox.com/s/2ibyxq98dfmnu1j/HI023_trimmed_2.fastq.gz?dl=0', 'HI023_trimmed_2.fastq.gz'], 'HI033_R1':['https://www.dropbox.com/s/58uia0fgwsucmds/HI033_trimmed_R1.fastq.gz?dl=0', 'HI033_trimmed_R1.fastq.gz'], 'HI033_R2':['https://www.dropbox.com/s/dujvkjphxypbb3o/HI033_trimmed_R2.fastq.gz?dl=0', 'HI033_trimmed_R2.fastq.gz'], 'mex986_R1':['https://www.dropbox.com/s/vi5x29okmhb8oi8/mex986_trimmed_R1.fastq.gz?dl=0', 'mex986_trimmed_R1.fastq.gz'], 'mex986_R2':['https://www.dropbox.com/s/n2eo6hxdvs2e7fr/mex986_trimmed_R2.fastq.gz?dl=0', 'mex986_trimmed_R2.fastq.gz'], 'mex919_R1':['https://www.dropbox.com/s/7oozv3zpljo13ob/mex919_trimmed_R1.fastq.gz?dl=0', 'mex919_trimmed_R1.fastq.gz'], 'mex919_R2':['https://www.dropbox.com/s/dqbs0ocx73kwqbh/mex919_trimmed_R2.fastq.gz?dl=0', 'mex919_trimmed_R2.fastq.gz'], 'mex915_R1':['https://www.dropbox.com/s/mzx2c49khxmt636/mex915_trimmed_R1.fastq.gz?dl=0', 'mex915_trimmed_R1.fastq.gz'], 'mex915_R2':['https://www.dropbox.com/s/n9pjb3rlql3mutk/mex915_trimmed_R2.fastq.gz?dl=0', 'mex915_trimmed_R2.fastq.gz'], 'mex536_R1':['https://www.dropbox.com/s/l5pfsodl6q23xfm/mex536_trimmed_R1.fastq.gz?dl=0', 'mex536_trimmed_R1.fastq.gz'], 'mex536_R2':['https://www.dropbox.com/s/wdq9nyjvd6sqs24/mex536_trimmed_R2.fastq.gz?dl=0', 'mex536_trimmed_R2.fastq.gz'], 'mex1527_R1':['https://www.dropbox.com/s/fj4eldilgas5cft/mex1527_trimmed_R1.fastq.gz?dl=0', 'mex1527_trimmed_R1.fastq.gz'], 'mex1527_R2':['https://www.dropbox.com/s/m1oqmq05bqli70j/mex1527_trimmed_R2.fastq.gz?dl=0', 'mex1527_trimmed_R2.fastq.gz']}
samples=['HH1','HH2','HH3','HH4','HH5','HH6','HH7','HH8','HH9','HH10','PL1','PL2','PL3','PL4','PL5','PL6','PL7','PL8','PL9','PL10','ESB1','ESB2','ESB3','ESB4','ESB5','ESB6','ESB7','ESB8','ESB9','ESB10','stm163','stm146','T9','T14','NJ203','NJ116','NJ1','HI023','HI033','mex986','mex919','mex915','mex536','mex1527']

link_dict_lb_2={'HH2_1_R1':['https://www.dropbox.com/s/kv07ovatdaennev/HH2_CAGATC_L005_R1_001.fastq.gz?dl=0','HH2_CAGATC_L005_R1_001.fastq.gz'],'HH2_1_R2':['https://www.dropbox.com/s/ocv9nmw8mrbxi6d/HH2_CAGATC_L005_R2_001.fastq.gz?dl=0','HH2_CAGATC_L005_R2_001.fastq.gz'],'HH2_2_R1':['https://www.dropbox.com/s/kcgzqk2kj7zxef2/HH2_CAGATC_L005_R1_002.fastq.gz?dl=0','HH2_CAGATC_L005_R1_002.fastq.gz'],'HH2_2_R2':['https://www.dropbox.com/s/6upkt817sm7eodb/HH2_CAGATC_L005_R2_002.fastq.gz?dl=0','HH2_CAGATC_L005_R2_002.fastq.gz'],'HH3_1_R1':['https://www.dropbox.com/s/17d1dc16kyxe0ph/HH3_ACTTGA_L005_R1_001.fastq.gz?dl=0','HH3_ACTTGA_L005_R1_001.fastq.gz'],'HH3_1_R2':['https://www.dropbox.com/s/jju72jws2wl2fsb/HH3_ACTTGA_L005_R2_001.fastq.gz?dl=0','HH3_ACTTGA_L005_R2_001.fastq.gz'],'HH3_2_R1':['https://www.dropbox.com/s/17d1dc16kyxe0ph/HH3_ACTTGA_L005_R1_001.fastq.gz?dl=0','HH3_ACTTGA_L005_R1_001.fastq.gz'],'HH3_2_R2':['https://www.dropbox.com/s/jju72jws2wl2fsb/HH3_ACTTGA_L005_R2_001.fastq.gz?dl=0','HH3_ACTTGA_L005_R2_001.fastq.gz'],'HH5_1_R1':['https://www.dropbox.com/s/v83oijg8roxp1hl/HH5_TAGCTT_L006_R1_001.fastq.gz?dl=0','HH5_TAGCTT_L006_R1_001.fastq.gz'],'HH5_1_R2':['https://www.dropbox.com/s/wc87ezu1vywk0zn/HH5_TAGCTT_L006_R2_001.fastq.gz?dl=0','HH5_TAGCTT_L006_R2_001.fastq.gz'],'HH5_2_R1':['https://www.dropbox.com/s/imm4loc6nwjbms5/HH5_TAGCTT_L006_R1_002.fastq.gz?dl=0','HH5_TAGCTT_L006_R1_002.fastq.gz'],'HH5_2_R2':['https://www.dropbox.com/s/0ie5b4u3tjrjfe3/HH5_TAGCTT_L006_R2_002.fastq.gz?dl=0','HH5_TAGCTT_L006_R2_002.fastq.gz'],'HH9_1_R1':['https://www.dropbox.com/s/1r3cxt9h15y5rv1/HH9_ATTCCT_L007_R1_001.fastq.gz?dl=0','HH9_ATTCCT_L007_R1_001.fastq.gz'],'HH9_1_R2':['https://www.dropbox.com/s/f4yorcjv4fe3yy6/HH9_ATTCCT_L007_R2_001.fastq.gz?dl=0','HH9_ATTCCT_L007_R2_001.fastq.gz'],'HH9_2_R1':['https://www.dropbox.com/s/9nslvjlalkerggh/HH9_ATTCCT_L007_R1_002.fastq.gz?dl=0','HH9_ATTCCT_L007_R1_002.fastq.gz'],'HH9_2_R2':['https://www.dropbox.com/s/4jhgtzk46hvf57w/HH9_ATTCCT_L007_R2_002.fastq.gz?dl=0','HH9_ATTCCT_L007_R2_002.fastq.gz'],'PL1_1_R1':['https://www.dropbox.com/s/14wxsl32don3ltd/PL1_ATCACG_L005_R1_001.fastq.gz?dl=0','PL1_ATCACG_L005_R1_001.fastq.gz'],'PL1_1_R2':['https://www.dropbox.com/s/utv5jn45c782qdb/PL1_ATCACG_L005_R2_001.fastq.gz?dl=0','PL1_ATCACG_L005_R2_001.fastq.gz'],'PL1_2_R1':['https://www.dropbox.com/s/0wa32dwtk8om008/PL1_ATCACG_L005_R1_002.fastq.gz?dl=0','PL1_ATCACG_L005_R1_002.fastq.gz'],'PL1_2_R2':['https://www.dropbox.com/s/3ma8ztsoiqn9gf2/PL1_ATCACG_L005_R2_002.fastq.gz?dl=0','PL1_ATCACG_L005_R2_002.fastq.gz'],'PL4_1_R1':['https://www.dropbox.com/s/7pypaknzivazzkv/PL4_TGACCA_L005_R1_001.fastq.gz?dl=0','PL4_TGACCA_L005_R1_001.fastq.gz'],'PL4_1_R2':['https://www.dropbox.com/s/6dwgx4a102xulpm/PL4_TGACCA_L005_R2_001.fastq.gz?dl=0','PL4_TGACCA_L005_R2_001.fastq.gz'],'PL4_2_R1':['https://www.dropbox.com/s/qv4fz6zjeoa4eob/PL4_TGACCA_L005_R1_002.fastq.gz?dl=0','PL4_TGACCA_L005_R1_002.fastq.gz'],'PL4_2_R2':['https://www.dropbox.com/s/bdxxnajcqdpssr8/PL4_TGACCA_L005_R2_002.fastq.gz?dl=0','PL4_TGACCA_L005_R2_002.fastq.gz'],'PL5_1_R1':['https://www.dropbox.com/s/9zrpr4bac1wc75q/PL5_ACAGTG_L006_R1_001.fastq.gz?dl=0','PL5_ACAGTG_L006_R1_001.fastq.gz'],'PL5_1_R2':['https://www.dropbox.com/s/sygh78hmjoki7ig/PL5_ACAGTG_L006_R2_001.fastq.gz?dl=0','PL5_ACAGTG_L006_R2_001.fastq.gz'],'PL5_2_R1':['https://www.dropbox.com/s/o6y0vyzpa7zs2o5/PL5_ACAGTG_L006_R1_002.fastq.gz?dl=0','PL5_ACAGTG_L006_R1_002.fastq.gz'],'PL5_2_R2':['https://www.dropbox.com/s/6oc2g0ia4tiw6b5/PL5_ACAGTG_L006_R2_002.fastq.gz?dl=0','PL5_ACAGTG_L006_R2_002.fastq.gz'],'PL6_1_R1':['https://www.dropbox.com/s/1phpfz4089b1vx9/PL6_CCGTCC_L006_R1_001.fastq.gz?dl=0','PL6_CCGTCC_L006_R1_001.fastq.gz'],'PL6_1_R2':['https://www.dropbox.com/s/ay62wei4cv1yeta/PL6_CCGTCC_L006_R2_001.fastq.gz?dl=0','PL6_CCGTCC_L006_R2_001.fastq.gz'],'PL6_2_R1':['https://www.dropbox.com/s/qhajeeivawa1nj9/PL6_CCGTCC_L006_R1_002.fastq.gz?dl=0','PL6_CCGTCC_L006_R1_002.fastq.gz'],'PL6_2_R2':['https://www.dropbox.com/s/55i6l9o992cennm/PL6_CCGTCC_L006_R2_002.fastq.gz?dl=0','PL6_CCGTCC_L006_R2_002.fastq.gz'],'PL7_1_R1':['https://www.dropbox.com/s/ssrtt2l9u09bvsc/PL7_GTCCGC_L006_R1_001.fastq.gz?dl=0','PL7_GTCCGC_L006_R1_001.fastq.gz'],'PL7_1_R2':['https://www.dropbox.com/s/ou2sp1q1c88yz5u/PL7_GTCCGC_L006_R2_001.fastq.gz?dl=0','PL7_GTCCGC_L006_R2_001.fastq.gz'],'PL7_2_R1':['https://www.dropbox.com/s/v4qr1ng8taluk4a/PL7_GTCCGC_L006_R1_002.fastq.gz?dl=0','PL7_GTCCGC_L006_R1_002.fastq.gz'],'PL7_2_R2':['https://www.dropbox.com/s/okbjpmcryprylmi/PL7_GTCCGC_L006_R2_002.fastq.gz?dl=0','PL7_GTCCGC_L006_R2_002.fastq.gz'],'PL9_1_R1':['https://www.dropbox.com/s/yshrv5ncpulm5o9/PL9_GTGGCC_L007_R1_001.fastq.gz?dl=0','PL9_GTGGCC_L007_R1_001.fastq.gz'],'PL9_1_R2':['https://www.dropbox.com/s/xnhfvf6b9x6x1p5/PL9_GTGGCC_L007_R1_002.fastq.gz?dl=0','PL9_GTGGCC_L007_R1_002.fastq.gz'],'PL9_2_R1':['https://www.dropbox.com/s/bvugi8kmmyrvgzo/PL9_GTGGCC_L007_R2_001.fastq.gz?dl=0','PL9_GTGGCC_L007_R2_001.fastq.gz'],'PL9_2_R2':['https://www.dropbox.com/s/qaqwb9a0pkd6qah/PL9_GTGGCC_L007_R2_002.fastq.gz?dl=0','PL9_GTGGCC_L007_R2_002.fastq.gz'],'PL10_1_R1':['https://www.dropbox.com/s/7b6hx5fgdvlk9np/PL10_GTTTCG_L007_R1_001.fastq.gz?dl=0','PL10_GTTTCG_L007_R1_001.fastq.gz'],'PL10_1_R2':['https://www.dropbox.com/s/zee5uyxbt94watk/PL10_GTTTCG_L007_R2_001.fastq.gz?dl=0','PL10_GTTTCG_L007_R2_001.fastq.gz'],'PL10_2_R1':['https://www.dropbox.com/s/t4qdtys93dxr30c/PL10_GTTTCG_L007_R1_002.fastq.gz?dl=0','PL10_GTTTCG_L007_R1_002.fastq.gz'],'PL10_2_R2':['https://www.dropbox.com/s/5azwidms8njr73y/PL10_GTTTCG_L007_R2_002.fastq.gz?dl=0','PL10_GTTTCG_L007_R2_002.fastq.gz'],'ESB1_1_R1':['https://www.dropbox.com/s/8lfodulfypu9qj8/ESB1_GGCTAC_L005_R1_001.fastq.gz?dl=0','ESB1_GGCTAC_L005_R1_001.fastq.gz'],'ESB1_1_R2':['https://www.dropbox.com/s/aibwjjxvdtsn0a5/ESB1_GGCTAC_L005_R2_001.fastq.gz?dl=0','ESB1_GGCTAC_L005_R2_001.fastq.gz'],'ESB1_2_R1':['https://www.dropbox.com/s/ag2l95m4n2iy1i4/ESB1_GGCTAC_L005_R1_002.fastq.gz?dl=0','ESB1_GGCTAC_L005_R1_002.fastq.gz'],'ESB1_2_R2':['https://www.dropbox.com/s/m53n3gp11sn1rum/ESB1_GGCTAC_L005_R2_002.fastq.gz?dl=0','ESB1_GGCTAC_L005_R2_002.fastq.gz'],'ESB2_1_R1':['https://www.dropbox.com/s/vj5eq67gdlrbmeb/ESB2_CTTGTA_L005_R1_001.fastq.gz?dl=0','ESB2_CTTGTA_L005_R1_001.fastq.gz'],'ESB2_1_R2':['https://www.dropbox.com/s/kh6ld8mdgqg94qb/ESB2_CTTGTA_L005_R1_002.fastq.gz?dl=0','ESB2_CTTGTA_L005_R1_002.fastq.gz'],'ESB2_2_R1':['https://www.dropbox.com/s/i7pjpcbwy94q0ox/ESB2_CTTGTA_L005_R2_001.fastq.gz?dl=0','ESB2_CTTGTA_L005_R2_001.fastq.gz'],'ESB2_2_R2':['https://www.dropbox.com/s/9wyj275m14kmafu/ESB2_CTTGTA_L005_R2_002.fastq.gz?dl=0','ESB2_CTTGTA_L005_R2_002.fastq.gz'],'ESB5_1_R1':['https://www.dropbox.com/s/2mab1h837wmy4ap/ESB5_ATGTCA_L006_R1_001.fastq.gz?dl=0','ESB5_ATGTCA_L006_R1_001.fastq.gz'],'ESB5_1_R2':['https://www.dropbox.com/s/l68wqlt3mdqtwcg/ESB5_ATGTCA_L006_R2_001.fastq.gz?dl=0','ESB5_ATGTCA_L006_R2_001.fastq.gz'],'ESB5_2_R1':['https://www.dropbox.com/s/fuk0xm43rye970u/ESB5_ATGTCA_L006_R1_002.fastq.gz?dl=0','ESB5_ATGTCA_L006_R1_002.fastq.gz'],'ESB5_2_R2':['https://www.dropbox.com/s/fkyxgxn4yu2d8q1/ESB5_ATGTCA_L006_R2_002.fastq.gz?dl=0','ESB5_ATGTCA_L006_R2_002.fastq.gz'],'ESB6_1_R1':['https://www.dropbox.com/s/1xnsdvvw5cbjgjb/ESB6_CTTGTA_L006_R1_001.fastq.gz?dl=0','ESB6_CTTGTA_L006_R1_001.fastq.gz'],'ESB6_1_R2':['https://www.dropbox.com/s/4h3bsfubdgxc3ak/ESB6_CTTGTA_L006_R2_001.fastq.gz?dl=0','ESB6_CTTGTA_L006_R2_001.fastq.gz'],'ESB6_2_R1':['https://www.dropbox.com/s/oreto7hy24om5ne/ESB6_CTTGTA_L006_R1_002.fastq.gz?dl=0','ESB6_CTTGTA_L006_R1_002.fastq.gz'],'ESB6_2_R2':['https://www.dropbox.com/s/mu32pqmlkq66el3/ESB6_CTTGTA_L006_R2_002.fastq.gz?dl=0','ESB6_CTTGTA_L006_R2_002.fastq.gz'],'ESB7_1_R1':['https://www.dropbox.com/s/xrxxjp79wb39kbo/ESB7_CTTGTA_L007_R1_001.fastq.gz?dl=0','ESB7_CTTGTA_L007_R1_001.fastq.gz'],'ESB7_1_R2':['https://www.dropbox.com/s/qereaub2xp255lv/ESB7_CTTGTA_L007_R2_001.fastq.gz?dl=0','ESB7_CTTGTA_L007_R2_001.fastq.gz'],'ESB7_2_R1':['https://www.dropbox.com/s/zqtyq58lpopmsar/ESB7_CTTGTA_L007_R1_002.fastq.gz?dl=0','ESB7_CTTGTA_L007_R1_002.fastq.gz'],'ESB7_2_R2':['https://www.dropbox.com/s/wkn1y40hrid3hec/ESB7_CTTGTA_L007_R2_002.fastq.gz?dl=0','ESB7_CTTGTA_L007_R2_002.fastq.gz']}
link_dict_lb_1={'HH1_R1':['https://www.dropbox.com/s/8btzg1mo9q7n1oo/HH1_GCCAAT_L005_R1_001.fastq.gz?dl=0', 'HH1_GCCAAT_L005_R1_001.fastq.gz'],'HH1_R2':['https://www.dropbox.com/s/p2p7f06k6tv9ylr/HH1_GCCAAT_L005_R2_001.fastq.gz?dl=0', 'HH1_GCCAAT_L005_R2_001.fastq.gz'],'HH4_R1':['https://www.dropbox.com/s/32uziv2uxvmshvl/HH4_GATCAG_L006_R1_001.fastq.gz?dl=0', 'HH4_GATCAG_L006_R1_001.fastq.gz'],'HH4_R2':['https://www.dropbox.com/s/18jeta1emw678ux/HH4_GATCAG_L006_R2_001.fastq.gz?dl=0', 'HH4_GATCAG_L006_R2_001.fastq.gz'],'HH6_R1':['https://www.dropbox.com/s/ke5qburf1w3ancv/HH6_CGTACG_L006_R1_001.fastq.gz?dl=0', 'HH6_CGTACG_L006_R1_001.fastq.gz'],'HH6_R2':['https://www.dropbox.com/s/u5przx6ylm8q89i/HH6_CGTACG_L006_R2_001.fastq.gz?dl=0', 'HH6_CGTACG_L006_R2_001.fastq.gz'],'HH7_R1':['https://www.dropbox.com/s/4z9l0zj13qvlsfk/HH7_GAGTGG_L006_R1_001.fastq.gz?dl=0', 'HH7_GAGTGG_L006_R1_001.fastq.gz'],'HH7_R2':['https://www.dropbox.com/s/qramrd46tzipu2f/HH7_GAGTGG_L006_R2_001.fastq.gz?dl=0', 'HH7_GAGTGG_L006_R2_001.fastq.gz'],'HH8_R1':['https://www.dropbox.com/s/v3qxuulaxyyf280/HH8_ACTGAT_L007_R1_001.fastq.gz?dl=0', 'HH8_ACTGAT_L007_R1_001.fastq.gz'],'HH8_R2':['https://www.dropbox.com/s/rcm5c0hkpsbkjq8/HH8_ACTGAT_L007_R2_001.fastq.gz?dl=0', 'HH8_ACTGAT_L007_R2_001.fastq.gz'],'HH10_R1':['https://www.dropbox.com/s/syeu04mljwmmk55/HH10_ATCACG_L007_R1_001.fastq.gz?dl=0', 'HH10_ATCACG_L007_R1_001.fastq.gz'],'HH10_R2':['https://www.dropbox.com/s/ab42503oqnzjotb/HH10_ATCACG_L007_R2_001.fastq.gz?dl=0', 'HH10_ATCACG_L007_R2_001.fastq.gz'],'PL2_R1':['https://www.dropbox.com/s/d3hx437d3qbtdmt/PL2_CGATGT_L005_R1_001.fastq.gz?dl=0', 'PL2_CGATGT_L005_R1_001.fastq.gz'],'PL2_R2':['https://www.dropbox.com/s/7cofey0oxbouy8e/PL2_CGATGT_L005_R2_001.fastq.gz?dl=0', 'PL2_CGATGT_L005_R2_001.fastq.gz'],'PL3_R1':['https://www.dropbox.com/s/dghx0l04xalvnio/PL3_TTAGGC_L005_R1_001.fastq.gz?dl=0', 'PL3_TTAGGC_L005_R1_001.fastq.gz'],'PL3_R2':['https://www.dropbox.com/s/7qnbrktwmsyju8a/PL3_TTAGGC_L005_R2_001.fastq.gz?dl=0', 'PL3_TTAGGC_L005_R2_001.fastq.gz'],'PL8_R1':['https://www.dropbox.com/s/59m4f1q12cqjiha/PL8_GTGAAA_L007_R1_001.fastq.gz?dl=0', 'PL8_GTGAAA_L007_R1_001.fastq.gz'] ,'PL8_R2':['https://www.dropbox.com/s/57r5s4zo61au7zs/PL8_GTGAAA_L007_R2_001.fastq.gz?dl=0', 'PL8_GTGAAA_L007_R2_001.fastq.gz'],'ESB3_R1':['https://www.dropbox.com/s/4ix2juq6yyu0431/ESB3_AGTCAA_L005_R1_001.fastq.gz?dl=0', 'ESB3_AGTCAA_L005_R1_001.fastq.gz'],'ESB3_R2':['https://www.dropbox.com/s/bd7axaqvh27hv2b/ESB3_AGTCAA_L005_R2_001.fastq.gz?dl=0', 'ESB3_AGTCAA_L005_R2_001.fastq.gz'],'ESB4_R1':['https://www.dropbox.com/s/e34un7x4he5amcb/ESB4_AGTTCC_L006_R1_001.fastq.gz?dl=0', 'ESB4_AGTTCC_L006_R1_001.fastq.gz'],'ESB4_R2':['https://www.dropbox.com/s/ppr7lzk26jy7mxp/ESB4_AGTTCC_L006_R2_001.fastq.gz?dl=0', 'ESB4_AGTTCC_L006_R2_001.fastq.gz'],'ESB8_R1':['https://www.dropbox.com/s/75qoo6e3wnicd2i/ESB8_TGACCA_L007_R1_001.fastq.gz?dl=0', 'ESB8_TGACCA_L007_R1_001.fastq.gz'],'ESB8_R2':['https://www.dropbox.com/s/4yr11r8zaikmgeg/ESB8_TGACCA_L007_R2_001.fastq.gz?dl=0', 'ESB8_TGACCA_L007_R2_001.fastq.gz'],'ESB9_R1':['https://www.dropbox.com/s/c1eyp0s2r0enwlp/ESB9_ACAGTG_L007_R1_001.fastq.gz?dl=0', 'ESB9_ACAGTG_L007_R1_001.fastq.gz'],'ESB9_R2':['https://www.dropbox.com/s/twf4bnc5hvm5j0u/ESB9_ACAGTG_L007_R2_001.fastq.gz?dl=0', 'ESB9_ACAGTG_L007_R2_001.fastq.gz'],'ESB10_R1':['https://www.dropbox.com/s/u6hy4gkjwa6u6c7/ESB10_GCCAAT_L007_R1_001.fastq.gz?dl=0', 'ESB10_GCCAAT_L007_R1_001.fastq.gz'] ,'ESB10_R2':['https://www.dropbox.com/s/uz3g35vwlcnjbho/ESB10_GCCAAT_L007_R2_001.fastq.gz?dl=0', 'ESB10_GCCAAT_L007_R2_001.fastq.gz'] ,'stm163_R1':['https://www.dropbox.com/s/2u525xxqpa64pgu/SRR1552223_1.fastq.gz?dl=0','stm163_R1.fastq.gz'],'stm163_R2':['https://www.dropbox.com/s/hbk8smo4zfrr28r/SRR1552223_2.fastq.gz?dl=0','stm163_R2.fastq.gz'],'stm146_R1':['https://www.dropbox.com/s/v9jr3nhtnszjqr9/SRR1552222_1.fastq.gz?dl=0','stm146_R1.fastq.gz'],'stm146_R2':['https://www.dropbox.com/s/w0xwqpo8j8pk0qf/SRR1552222_2.fastq.gz?dl=0','stm146_R2.fastq.gz'],'T9_R1':['https://www.dropbox.com/s/enk0hv2oyrwmrej/SRR1549529_1.fastq.gz?dl=0','T9_R1.fastq.gz'] ,'T9_R2':['https://www.dropbox.com/s/eyhmsor6cu4hsc4/SRR1549529_2.fastq.gz?dl=0','T9_R2.fastq.gz'],'T14_R1':['https://www.dropbox.com/s/4ykzq1ifqi3m98d/SRR1549528_1.fastq.gz?dl=0','T14_R1.fastq.gz'],'T14_R2':['https://www.dropbox.com/s/dbpjf288bsovxp1/SRR1549528_2.fastq.gz?dl=0','T14_R2.fastq.gz'],'NJ203_R1':['https://www.dropbox.com/s/23zbqsnd85salto/SRR1548574_1.fastq.gz?dl=0','NJ203_R1.fastq.gz'],'NJ203_R2':['https://www.dropbox.com/s/aqpnbnjzew9gird/SRR1548574_2.fastq.gz?dl=0','NJ203_R2.fastq.gz'],'NJ116_R1':['https://www.dropbox.com/s/yq643jcnsc3mz6t/SRR1548573_1.fastq.gz?dl=0', 'NJ116_R1.fastq.gz'],'NJ116_R2':['https://www.dropbox.com/s/td37rxh9a3jhesb/SRR1548573_2.fastq.gz?dl=0', 'NJ116_R2.fastq.gz'],'NJ1_R1':['https://www.dropbox.com/s/aehx9gikdm8fuh2/SRR1548572_1.fastq.gz?dl=0','NJ1_R1.fastq.gz'],'NJ1_R2':['https://www.dropbox.com/s/hc7qy4lwa5paukb/SRR1548572_2.fastq.gz?dl=0','NJ1_R2.fastq.gz'],'HI023_R1':['https://www.dropbox.com/s/ygbvczsrq88vjxm/SRR1548571_1.fastq.gz?dl=0','HI023_R1.fastq.gz'],'HI023_R2':['https://www.dropbox.com/s/8j1hpn8bxebmzdn/SRR1548571_2.fastq.gz?dl=0','HI023_R2.fastq.gz'],'HI033_R1':['https://www.dropbox.com/s/c71ey1yf46c0ld5/SRR1548576_1.fastq.gz?dl=0','HI033_R1.fastq.gz'],'HI033_R2':['https://www.dropbox.com/s/6h7qjopwq93pq7e/SRR1548576_2.fastq.gz?dl=0','HI033_R2.fastq.gz'],'mex986_R1':['https://www.dropbox.com/s/lsznzg1jc3iurdq/SRR1552209_1.fastq.gz?dl=0','mex986_R1.fastq.gz'],'mex986_R2':['https://www.dropbox.com/s/kf6dvvlrzut2rc5/SRR1552209_2.fastq.gz?dl=0','mex986_R2.fastq.gz'],'mex919_R1':['https://www.dropbox.com/s/b47c6l4egqt9e4d/SRR1552208_1.fastq.gz?dl=0', 'mex919_R1.fastq.gz'],'mex919_R2':['https://www.dropbox.com/s/j5s63oaqu4enq5y/SRR1552208_2.fastq.gz?dl=0', 'mex919_R2.fastq.gz'],'mex915_R1':['https://www.dropbox.com/s/4ym0g0fjpd5w23z/SRR1552207_1.fastq.gz?dl=0','mex915_R1.fastq.gz'],'mex915_R2':['https://www.dropbox.com/s/jfv31ptntmj0xfl/SRR1552207_2.fastq.gz?dl=0','mex915_R2.fastq.gz'],'mex536_R1':['https://www.dropbox.com/s/yzpi8yj9zj1dbqu/SRR1552206_1.fastq.gz?dl=0','mex536_R1.fastq.gz'],'mex536_R2':['https://www.dropbox.com/s/u8m8owi776e6eg3/SRR1552206_2.fastq.gz?dl=0','mex536_R2.fastq.gz'],'mex1527_R1':['https://www.dropbox.com/s/tj1gxjw8kq63zrh/SRR1552204_1.fastq.gz?dl=0','mex1527_R1.fastq.gz'],'mex1527_R2':['https://www.dropbox.com/s/wls3afufglafter/SRR1552204_2.fastq.gz?dl=0','mex1527_R2.fastq.gz']}


for i in samples_1_libraries:
	fastq_1 = link_dict_lb_1[i+'_R1']
	fastq_2 = link_dict_lb_1[i+'_R2']
	Sample_ID=str(i)
	out_dir='/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/'
	genome='/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa'
	outfile='/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/'+Sample_ID+'.sh'
	run_BWA_GATK_lib_1(Sample_ID, genome, out_dir, fastq_1, fastq_2, outfile)


for i in samples_2_libraries:
	fastq_1 = link_dict_lb_2[i+'_1_R1']
	fastq_2 = link_dict_lb_2[i+'_1_R2']
	fastq_3 = link_dict_lb_2[i+'_2_R1']
	fastq_4 = link_dict_lb_2[i+'_2_R2']	
	Sample_ID=str(i)
	out_dir='/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/'
	genome='/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa'
	outfile='/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/'+Sample_ID+'.sh'
	run_BWA_GATK_lib_2(Sample_ID, genome, out_dir, fastq_1, fastq_2, fastq_3, fastq_4, outfile)



run_BCF_tools_snp_calling('stm146', '/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa', '/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/', '/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/'+'stm146_Bcftools'+'.sh')
run_BCF_tools_snp_calling('HH9', '/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa', '/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/', '/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/'+'HH9_Bcftools'+'.sh')


Knownsites="golden_set.vcf"
run_final_bam_and_call_gVCF(samples, Knownsites)


