# Here we provide the scripts to calculate the nucleotide diversity of  Eastern and Western monarchs for different genomic positions. 

Genomic positions were classified as Intergenic, Intronic, 1st, 2nd, 3rd codon positions and 4-fold degenerate sites (4D).


Before you use this code to callculate nucleotide diversity you will need the follwing files. 
1) a BED file with the genomic windows for which you would like to calculate the nucleotide diversity. (CHROM, BIN_START, BIN_END). (Windows.bed)
2) genome wide coverage of all the samples beeing used. it can be obtained using the below code. 
```
genomeCoverageBed -d -ibam mex536.final.bam  > mex536.'genomeCoverage.txt' &
genomeCoverageBed -d -ibam mex919.final.bam  > mex919.'genomeCoverage.txt' &
genomeCoverageBed -d -ibam mex986.final.bam  > mex986.'genomeCoverage.txt' &
genomeCoverageBed -d -ibam NJ116.final.bam  > NJ116.'genomeCoverage.txt' &
genomeCoverageBed -d -ibam PL2.final.bam  > PL2.'genomeCoverage.txt' &
genomeCoverageBed -d -ibam PL4.final.bam  > PL4.'genomeCoverage.txt' &
genomeCoverageBed -d -ibam PL5.final.bam  > PL5.'genomeCoverage.txt' &
genomeCoverageBed -d -ibam stm146.final.bam  > stm146.'genomeCoverage.txt' &
genomeCoverageBed -d -ibam T14.final.bam  > T14.'genomeCoverage.txt' &
genomeCoverageBed -d -ibam T9.final.bam  > T9.'genomeCoverage.txt' &
```
3)Depth of all SNPs/genotypes obtained using the below command. (.gdepth file)
```
vcftools --gzvcf recal_snps_all_20200810.vcf.gz --keep ./lists/randomly_selected/East.txt --remove-filtered-all --geno-depth --out East &
vcftools --gzvcf recal_snps_all_20200810.vcf.gz --keep ./lists/randomly_selected/West.txt --remove-filtered-all --geno-depth --out West &
```

commands to calculate the nucleotide diversity.
```
#j is coverage
#i is the subdivisions in the genome
for j in 2 3 4
do
for i in {0..19}
do 
echo 'Section of the genome: '$i  'Coverage: '$j
python get_N_sites_V1kb.py -n 8 -i East_bedfiles -o  'East_10000_'$j'X'  -l 10000 -p $i -c $j -R '/work/smalab/Venkat/genome/Danaus_plexippus_v3_-_scaffolds.fa' -t 'all' &
python get_N_sites_V1kb.py -n 8 -i West_bedfiles -o  'West_10000_'$j'X'  -l 10000 -p $i -c $j -R '/work/smalab/Venkat/genome/Danaus_plexippus_v3_-_scaffolds.fa' -t 'all' &
sleep 10m
done
done

for j in 2 3
do
for i in {0..19}
do 
echo $i $j
python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/East.frq' -S 'East_10000_'$j'X_all_'$i'_N_sites' -I 'all' -P '/scratch/vt20265/bam_files' -Sec $i -D '/scratch/vt20265/East.gdepth' -c $j  -Ind "/scratch/vt20265/lists/randomly_selected/East.txt" &
python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/West.frq' -S 'West_10000_'$j'X_all_'$i'_N_sites' -I 'all' -P '/scratch/vt20265/bam_files' -Sec $i -D '/scratch/vt20265/West.gdepth' -c $j  -Ind "/scratch/vt20265/lists/randomly_selected/West.txt" &
sleep 5m
done
done
```

# Monarch_Pi_all_sites.py
The script Monarch_Pi_all_sites.py generates the data base of nucleotide diversity seperatly for Intergenic, Intronic, 1st, 2nd, 3rd codon positions and 4-fold degenerate sites (4D). These results are provided in Pi_data.csv

# Monarch_plots_pi.py
This script will give the folowing plot
![alt text](https://github.com/venta380/Monarch_genomics/blob/master/Nucleotide%20diversity/Screenshot%202020-01-27%2012.34.54.png "Logo Title Text 1")

# Monarch_pi_tw_td_sfs_all_sites_1977.py
This script will calculate the nucleotide diversity (π), Watterson’s θ and Tajima’s D for the specific sites (All, intergenic, codon positions, 0fold and 4fold) in the genome. 

