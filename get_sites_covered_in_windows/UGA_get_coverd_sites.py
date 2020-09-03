import sys
import string
import numpy as np
import fileinput
import pandas
import personal_popgen
import argparse
import os
import functools

pwd='/scratch/vt20265/bam_files/'
os.chdir(pwd)


def get_args():
    parser = argparse.ArgumentParser(description='''outputs sites covered 1X in all Induveduals in population pairs

        -i file with list of paths to population lists
        -o output_prefix in database
        -p which part of the scafolds scaf1 scaf2 scaf3
        -n number of induvedluas in the set
        -l length of the window
        -t type of sites ('intergenic'/'all')
        -R refrence fasta file
        ''')

    parser.add_argument('-i', '--input', type=str, required=True)
    parser.add_argument('-o', '--output', type=str, required=True)
    parser.add_argument('-p', '--part', type=str, required=True)
    parser.add_argument('-n', '--num', type=int, required=True)
    parser.add_argument('-l', '--length', type=int, required=True)
    parser.add_argument('-t', '--type', type=str, required=True)
    parser.add_argument('-r', '--Ref', type=str, required=True)

    return parser.parse_args()

args = get_args()
coverage_list_pop1=args.input
coverage_list_pop1=[stuff.strip() for stuff in open(coverage_list_pop1)]
output=args.output
inds=int(args.num)
part=int(args.part)
type_of_sites=args.type
refrence_file=args.Ref

#coverage_list_pop1='1977_bedfiles'
#coverage_list_pop1=[stuff.strip() for stuff in open(coverage_list_pop1)]
#output='1977_10000'
#inds=8
#part=0
#type_of_sites='intergenic'
#refrence_file='/work/smalab/Venkat/genome/Danaus_plexippus_v3_-_scaffolds.fa'
############################# making windows ################################
#chromosome=pandas.read_table("/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/Fst_window_1000.windowed.weir.fst")[['CHROM','BIN_START','BIN_END']]


#chromosome=personal_popgen.make_windows('/work/smalab/Venkat/genome/Danaus_plexippus_v3_-_scaffolds.fa',args.length)
chromosome=pandas.read_csv("/scratch/vt20265/bam_files/windows.bed", sep='\t', names=['CHROM','BIN_START','BIN_END'], header=None)

temp=[i[0] for i in zip(np.array_split(list(set(chromosome.CHROM)),20)[int(part)])]
chromosome=chromosome[chromosome['CHROM'].isin(temp)]



#############################################################################



def  assess_coverage_filter_7X_induv(list_1, inds, scaffolds):
        headera=['CHROM', 'POS', 'COV0']
        keys=['COV0']
        list_all=[]
        for i in range(0,len(list_1)):
                if i == 0:
                        n=pandas.read_table(list_1[i], header=None, engine='c', sep='\t')
                        n.columns = headera
                        n=n[n['CHROM'].isin(scaffolds)]
                        list_all=[n]
                else:
                        n=pandas.read_table(list_1[i], header=None, engine='c', sep='\t')
                        headera=['CHROM', 'POS', 'COV'+str(i)]
                        n.columns = headera
                        n=n[n['CHROM'].isin(scaffolds)]
                        list_all.append(n)
                        keys.append('COV'+str(i))
                        nextone=personal_popgen.join_bam_coverage(list_all)
                        del list_all
                        list_all=[nextone]
        del list_all
        temp=nextone[keys]
        nextone=nextone[temp[temp >= 1].count(axis=1) >= inds]
        del temp
        return nextone[['CHROM','POS']]

stats_dir="/scratch/vt20265/bam_files/"
list_1=[]
list_1.append(assess_coverage_filter_7X_induv(coverage_list_pop1, inds, temp))

print ("loaded all files")
    
if type_of_sites == 'all'
    df_pass=functools.reduce(lambda x, y: pandas.merge(x, y, on=['CHROM','POS'], how='inner'), list_1)
    print ("merged all")
    del list_1
    df_pass['BIN_START']=(np.floor(df_pass['POS']/args.length)*args.length)+1
    df_pass['BIN_END']=df_pass['BIN_START']+(args.length-1)
    df_merge = pandas.merge(df_pass, chromosome, on=['CHROM','BIN_START','BIN_END'],how='inner')
    df_merge = df_merge.query('BIN_START <= POS and BIN_END >= POS')
    print ("computed windows")
    del df_pass
    out=pandas.DataFrame({output+'_N_sites' :df_merge.groupby(['CHROM','BIN_START','BIN_END']).size()}).reset_index()
    print ("computed output")
    out.to_csv('/scratch/vt20265/bam_files/'+output+'_'+str(part)+'_N_sites')
    print ("completed")
elif type_of_sites == 'intergenic'
    fasta_DPlex=personal_popgen.fasta_dict(refrence_file)
    codon_df_1=pandas.read_csv('codon_df_1.csv',names=['CHROM','POS'], sep=' ',header=None)
    codon_df_1=codon_df_1[codon_df_1['CHROM'].isin(temp)]
    codon_df_2=pandas.read_csv('codon_df_2.csv',names=['CHROM','POS'], sep=' ',header=None)
    codon_df_2=codon_df_2[codon_df_2['CHROM'].isin(temp)]
    codon_df_3=pandas.read_csv('codon_df_3.csv',names=['CHROM','POS'], sep=' ',header=None)
    codon_df_3=codon_df_3[codon_df_3['CHROM'].isin(temp)]
    introns=pandas.read_csv('introns.csv',names=['CHROM','POS'], sep=' ',header=None)
    introns=introns[introns['CHROM'].isin(temp)]
    intergenic_exclude=pandas.concat([codon_df_1[['CHROM','POS']], codon_df_2[['CHROM','POS']], codon_df_3[['CHROM','POS']], introns[['CHROM','POS']]], ignore_index=True)
    intergenic_exclude=intergenic_exclude[['CHROM','POS']]
    intergenic=pandas.DataFrame(columns=['CHROM','POS'])
    for i in temp:
        chrom_len=len(fasta_DPlex[i])
        all_site=list(range(1,chrom_len+1))
        exlude_sites=list(intergenic_exclude[(intergenic_exclude.CHROM == i )]['POS'])
        temp=pandas.DataFrame({'CHROM':[i]*len(all_site), 'POS':all_site})
        temp=temp[~temp['POS'].isin(exlude_sites)]
        intergenic=intergenic.append(temp, ignore_index=True)
    all_sites=intergenic
    df_pass=functools.reduce(lambda x, y: pandas.merge(x, y, on=['CHROM','POS'], how='inner'), list_1)
    print ("merged all")
    df_pass=df_pass.merge(all_sites, on=['CHROM','POS'], how='inner')
    del list_1
    df_pass['BIN_START']=(np.floor(df_pass['POS']/args.length)*args.length)+1
    df_pass['BIN_END']=df_pass['BIN_START']+(args.length-1)
    df_merge = pandas.merge(df_pass, chromosome, on=['CHROM','BIN_START','BIN_END'],how='inner')
    df_merge = df_merge.query('BIN_START <= POS and BIN_END >= POS')
    print ("computed windows")
    del df_pass
    out=pandas.DataFrame({output+'_N_sites' :df_merge.groupby(['CHROM','BIN_START','BIN_END']).size()}).reset_index()
    print ("computed output")
    out.to_csv('/scratch/vt20265/bam_files/'+output+'_intergenic_'+str(part)+'_N_sites')
    print ("completed")




