#PBS -S /bin/bash
#PBS -q batch
#PBS -l nodes=1:ppn=10
#PBS -N Filter_COV
#PBS -l walltime=60:00:00
#PBS -o Filter_COV.out
#PBS -e Filter_COV.err
#PBS -l mem=100gb



cd /scratch/vt20265/bam_files/


#genomeCoverageBed -d -ibam 02_77.final.bam  > 02_77.'genomeCoverage.txt' &
#genomeCoverageBed -d -ibam 03_77.final.bam  > 03_77.'genomeCoverage.txt' &
#genomeCoverageBed -d -ibam 04_77.final.bam  > 04_77.'genomeCoverage.txt' &
#genomeCoverageBed -d -ibam 05_77.final.bam  > 05_77.'genomeCoverage.txt' &
#genomeCoverageBed -d -ibam 07_77.final.bam  > 07_77.'genomeCoverage.txt' &
#genomeCoverageBed -d -ibam 08_77.final.bam  > 08_77.'genomeCoverage.txt' &
#genomeCoverageBed -d -ibam 10_77.final.bam  > 10_77.'genomeCoverage.txt' &
#genomeCoverageBed -d -ibam 12_77.final.bam  > 12_77.'genomeCoverage.txt' &
#genomeCoverageBed -d -ibam ESB10.final.bam  > ESB10.'genomeCoverage.txt' &
#genomeCoverageBed -d -ibam ESB6.final.bam  > ESB6.'genomeCoverage.txt' &
#wait
#genomeCoverageBed -d -ibam ESB9.final.bam  > ESB9.'genomeCoverage.txt' &
#genomeCoverageBed -d -ibam HH2.final.bam  > HH2.'genomeCoverage.txt' &
#genomeCoverageBed -d -ibam HH5.final.bam  > HH5.'genomeCoverage.txt' &
#genomeCoverageBed -d -ibam HI023.final.bam  > HI023.'genomeCoverage.txt' &
#genomeCoverageBed -d -ibam mex536.final.bam  > mex536.'genomeCoverage.txt' &
#genomeCoverageBed -d -ibam mex919.final.bam  > mex919.'genomeCoverage.txt' &
#genomeCoverageBed -d -ibam mex986.final.bam  > mex986.'genomeCoverage.txt' &
#genomeCoverageBed -d -ibam NJ116.final.bam  > NJ116.'genomeCoverage.txt' &
#genomeCoverageBed -d -ibam PL2.final.bam  > PL2.'genomeCoverage.txt' &
#genomeCoverageBed -d -ibam PL4.final.bam  > PL4.'genomeCoverage.txt' &
#genomeCoverageBed -d -ibam PL5.final.bam  > PL5.'genomeCoverage.txt' &
#genomeCoverageBed -d -ibam stm146.final.bam  > stm146.'genomeCoverage.txt' &
#genomeCoverageBed -d -ibam T14.final.bam  > T14.'genomeCoverage.txt' &
#genomeCoverageBed -d -ibam T9.final.bam  > T9.'genomeCoverage.txt' &
#wait



#echo "./lists/East" >> East_West
#echo "./lists/West" >> East_West

#West_bedfiles
#East_bedfiles
#1977_bedfiles

#for i in {0..19}; 
#do echo $i; 
#python get_N_sites_V1kb.py -n 8 -i 1977_bedfiles -o  1977_10000  -l 10000 -p $i -R '/work/smalab/Venkat/genome/Danaus_plexippus_v3_-_scaffolds.fa' -t 'all' & 
#sleep 10m
#python get_N_sites_V1kb.py -n 8 -i West_bedfiles -o  West_10000  -l 10000 -p $i -R '/work/smalab/Venkat/genome/Danaus_plexippus_v3_-_scaffolds.fa' -t 'all' & 
#sleep 10m
#python get_N_sites_V1kb.py -n 8 -i East_bedfiles -o  East_10000  -l 10000 -p $i -R '/work/smalab/Venkat/genome/Danaus_plexippus_v3_-_scaffolds.fa' -t 'all' & 
#sleep 10m
#done
#
##python merge_subprocesses.py
#
#for i in {0..19}; 
#do echo $i; 
#python get_N_sites_V1kb.py -n 8 -i 1977_bedfiles -o  1977_10000  -l 10000 -p $i -R '/work/smalab/Venkat/genome/Danaus_plexippus_v3_-_scaffolds.fa' -t 'intergenic' & 
#sleep 10m
#python get_N_sites_V1kb.py -n 8 -i West_bedfiles -o  West_10000  -l 10000 -p $i -R '/work/smalab/Venkat/genome/Danaus_plexippus_v3_-_scaffolds.fa' -t 'intergenic' & 
#sleep 10m
#python get_N_sites_V1kb.py -n 8 -i East_bedfiles -o  East_10000  -l 10000 -p $i -R '/work/smalab/Venkat/genome/Danaus_plexippus_v3_-_scaffolds.fa' -t 'intergenic' & 
#sleep 10m
#done

#python merge_subprocesses.py



#python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/1977.frq' -S 'codon_df_1.csv' -I 'codon1' -P '/scratch/vt20265/bam_files' -Sec $i
#python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/1977.frq' -S 'codon_df_2.csv' -I 'codon2' -P '/scratch/vt20265/bam_files' -Sec $i
#python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/1977.frq' -S 'codon_df_3.csv' -I 'codon3' -P '/scratch/vt20265/bam_files' -Sec $i
#python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/1977.frq' -S 'codon_df_4D.csv' -I '4fold' -P '/scratch/vt20265/bam_files' -Sec $i
#python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/1977.frq' -S 'introns.csv' -I 'introns' -P '/scratch/vt20265/bam_files' -Sec $i
#python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/1977.frq' -S 'intergenic' -I 'intergenic' -P '/scratch/vt20265/bam_files' -Sec $i
#
#python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/West.frq' -S 'codon_df_1.csv' -I 'codon1' -P '/scratch/vt20265/bam_files' -Sec $i
#python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/West.frq' -S 'codon_df_2.csv' -I 'codon2' -P '/scratch/vt20265/bam_files' -Sec $i
#python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/West.frq' -S 'codon_df_3.csv' -I 'codon3' -P '/scratch/vt20265/bam_files' -Sec $i
#python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/West.frq' -S 'codon_df_4D.csv' -I '4fold' -P '/scratch/vt20265/bam_files' -Sec $i
#python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/West.frq' -S 'introns.csv' -I 'introns' -P '/scratch/vt20265/bam_files' -Sec $i
#python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/West.frq' -S 'intergenic' -I 'intergenic' -P '/scratch/vt20265/bam_files' -Sec $i
#
#python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/East.frq' -S 'codon_df_1.csv' -I 'codon1' -P '/scratch/vt20265/bam_files' -Sec $i
#python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/East.frq' -S 'codon_df_2.csv' -I 'codon2' -P '/scratch/vt20265/bam_files' -Sec $i
#python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/East.frq' -S 'codon_df_3.csv' -I 'codon3' -P '/scratch/vt20265/bam_files' -Sec $i
#python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/East.frq' -S 'codon_df_4D.csv' -I '4fold' -P '/scratch/vt20265/bam_files' -Sec $i
#python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/East.frq' -S 'introns.csv' -I 'introns' -P '/scratch/vt20265/bam_files' -Sec $i
#python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/East.frq' -S 'intergenic' -I 'intergenic' -P '/scratch/vt20265/bam_files' -Sec $i


for i in {0..19}; 
do echo $i; 
python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/East.frq' -S 'East_10000_all_'$i'_N_sites' -I 'all' -P '/scratch/vt20265/bam_files' -Sec $i &
python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/West.frq' -S 'West_10000_all_'$i'_N_sites' -I 'all' -P '/scratch/vt20265/bam_files' -Sec $i &
python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/1977.frq' -S '1977_10000_all_'$i'_N_sites' -I 'all' -P '/scratch/vt20265/bam_files' -Sec $i &
sleep 10m
python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/East.frq' -S 'East_10000_intergenic_'$i'_N_sites' -I 'intergenic' -P '/scratch/vt20265/bam_files' -Sec $i &
python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/West.frq' -S 'West_10000_intergenic_'$i'_N_sites' -I 'intergenic' -P '/scratch/vt20265/bam_files' -Sec $i &
python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/1977.frq' -S '1977_10000_intergenic_'$i'_N_sites' -I 'intergenic' -P '/scratch/vt20265/bam_files' -Sec $i &
sleep 10m
python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/East.frq' -S 'introns.csv' -I 'introns' -P '/scratch/vt20265/bam_files' -Sec $i &
python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/West.frq' -S 'introns.csv' -I 'introns' -P '/scratch/vt20265/bam_files' -Sec $i &
python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/1977.frq' -S 'introns.csv' -I 'introns' -P '/scratch/vt20265/bam_files' -Sec $i &
sleep 10m
python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/East.frq' -S 'codon_df_1.csv' -I 'codon1' -P '/scratch/vt20265/bam_files' -Sec $i &
python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/West.frq' -S 'codon_df_1.csv' -I 'codon1' -P '/scratch/vt20265/bam_files' -Sec $i &
python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/1977.frq' -S 'codon_df_1.csv' -I 'codon1' -P '/scratch/vt20265/bam_files' -Sec $i &
sleep 10m
python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/East.frq' -S 'codon_df_2.csv' -I 'codon2' -P '/scratch/vt20265/bam_files' -Sec $i &
python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/West.frq' -S 'codon_df_2.csv' -I 'codon2' -P '/scratch/vt20265/bam_files' -Sec $i &
python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/1977.frq' -S 'codon_df_2.csv' -I 'codon2' -P '/scratch/vt20265/bam_files' -Sec $i &
sleep 10m
python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/East.frq' -S 'codon_df_3.csv' -I 'codon3' -P '/scratch/vt20265/bam_files' -Sec $i &
python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/West.frq' -S 'codon_df_3.csv' -I 'codon3' -P '/scratch/vt20265/bam_files' -Sec $i &
python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/1977.frq' -S 'codon_df_3.csv' -I 'codon3' -P '/scratch/vt20265/bam_files' -Sec $i &
sleep 10m
python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/East.frq' -S 'codon_df_4D.csv' -I '4fold' -P '/scratch/vt20265/bam_files' -Sec $i &
python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/West.frq' -S 'codon_df_4D.csv' -I '4fold' -P '/scratch/vt20265/bam_files' -Sec $i &
python Monarch_pi_tw_td_sfs_all_sites_1977.py -F '/scratch/vt20265/1977.frq' -S 'codon_df_4D.csv' -I '4fold' -P '/scratch/vt20265/bam_files' -Sec $i &
sleep 10m
done


