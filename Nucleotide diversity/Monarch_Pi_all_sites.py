import sys
import os
import string
import pandas
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
import math
import time
import sys
import personal_popgen
from itertools import izip



sns.set(font_scale=1.5)
sns.set_style("whitegrid", {'axes.grid' : False})

pwd='/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/'
os.chdir(pwd)




fasta=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa')

chromosome = pandas.DataFrame(columns=['CHROM','BIN_START','BIN_END'])



for i in fasta.keys():
        temp={}
        if len(fasta[i]) >= 10000:
                list_start=range(1,len(fasta[i]),10000)
                list_start.pop()
                list_end=range(10000,len(fasta[i]),10000)
        else:
                list_start=[1]
                list_end=[10000]
        chrom=[i]*len(list_end)
        temp['CHROM']=pandas.Series(chrom)
        temp['BIN_START']=pandas.Series(list_start)
        temp['BIN_END']=pandas.Series(list_end)
        temp_df = pandas.DataFrame(temp)
        chromosome=personal_popgen.join_data_base([chromosome,temp_df])


temp=[i[0] for i in izip(np.array_split(list(set(chromosome.CHROM)),10)[int(sys.argv[1])])]

table_windows_2=chromosome[chromosome['CHROM'].isin(temp)]


def join_data_bases(lists):
    for i in range(1,len(lists)):
        j=i-1
        if i == 1:
            new_2=pandas.merge(lists[j], lists[i], on=['CHROM','BIN_START','BIN_END','POS'], how='outer')
        else:
            new_2=pandas.merge(nextone, lists[i], on=['CHROM','BIN_START','BIN_END','POS'], how='outer')
        nextone=new_2
    return nextone


def get_site_freq_spectrum(goods,codon_df_1):
    codon_df_1_temp=pandas.merge(goods,codon_df_1, on=['CHROM','POS'], how='inner')
    codon_df_1_temp.minor_freq=codon_df_1_temp.minor_freq.round(2)
    spectrum=codon_df_1_temp.groupby('minor_freq')['minor_freq'].count()
    return spectrum


codon_df_0_exclude=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/codon_df_0_exclude.csv.gz',names=['CHROM','POS'], sep=' ',header=None)
introns=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/introns.csv.gz',names=['POS','CHROM'], sep=' ',header=None)
codon_df_1=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/codon_df_1.csv.gz',names=['CHROM','POS','major_allele'], sep=' ',header=None)
codon_df_2=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/codon_df_2.csv.gz',names=['CHROM','POS','major_allele'], sep=' ',header=None)
codon_df_3=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/codon_df_3.csv.gz',names=['CHROM','POS','major_allele'], sep=' ',header=None)
codon_df_4d=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/codon_df_4d.csv.gz',names=['CHROM','POS'], sep=' ',header=None)

codon_df_0_exclude=codon_df_0_exclude[codon_df_0_exclude['CHROM'].isin(temp)]
introns=introns[introns['CHROM'].isin(temp)]
codon_df_1=codon_df_1[codon_df_1['CHROM'].isin(temp)]
codon_df_2=codon_df_2[codon_df_2['CHROM'].isin(temp)]
codon_df_3=codon_df_3[codon_df_3['CHROM'].isin(temp)]
codon_df_4d=codon_df_4d[codon_df_4d['CHROM'].isin(temp)]


codon_df_0d=pandas.concat([codon_df_1, codon_df_2], ignore_index=True)
codon_df_0d=codon_df_0d.merge(codon_df_0_exclude, how='outer', indicator=True)
codon_df_0d=codon_df_0d.query('_merge == "left_only"').drop('_merge', 1).drop('major_allele', 1)

intergenic_final=pandas.DataFrame(columns=['CHROM','POS'])


table_windows_2['reference']=0.0
table_windows_2['intron']=0.0
table_windows_2['codon_1']=0.0
table_windows_2['codon_2']=0.0
table_windows_2['codon_3']=0.0
table_windows_2['codon_4D']=0.0
table_windows_2['intergenic']=0.0


for window_index, window in table_windows_2.iterrows():
    fraction=len(str(fasta[window.CHROM][int(window.BIN_START-1):int(window.BIN_END-1)]))
    sequnce=np.array(fasta[window.CHROM])
    introns_seq_1=introns[(introns['CHROM']==window.CHROM) & (introns['POS']>=window.BIN_START) & (introns['POS']<=window.BIN_END) ]
    introns_seq=list(introns_seq_1['POS']-1)
    introns_seq=sequnce[introns_seq]
    introns_seq=float(len(''.join(introns_seq)))
    CD_1_seq_1=codon_df_1[(codon_df_1['CHROM']==window.CHROM) & (codon_df_1['POS']>=window.BIN_START) & (codon_df_1['POS']<=window.BIN_END) ]
    CD_1_seq=list(CD_1_seq_1['POS']-1)
    if len(CD_1_seq) != 0 and max(CD_1_seq) <= len(sequnce):
        CD_1_seq=sequnce[CD_1_seq]
        CD_1_seq=float(len(''.join(CD_1_seq)))
    else:
        CD_1_seq=np.nan
    CD_2_seq_1=codon_df_2[(codon_df_2['CHROM']==window.CHROM) & (codon_df_2['POS']>=window.BIN_START) & (codon_df_2['POS']<=window.BIN_END) ]
    CD_2_seq=list(CD_2_seq_1['POS']-1)
    if len(CD_2_seq) != 0 and max(CD_2_seq) <= len(sequnce):
        CD_2_seq=sequnce[CD_2_seq]
        CD_2_seq=float(len(''.join(CD_2_seq)))
    else:
        CD_2_seq=np.nan
    CD_3_seq_1=codon_df_3[(codon_df_3['CHROM']==window.CHROM) & (codon_df_3['POS']>=window.BIN_START) & (codon_df_3['POS']<=window.BIN_END) ]
    CD_3_seq=list(CD_3_seq_1['POS']-1)
    if len(CD_3_seq) != 0 and max(CD_3_seq) < len(sequnce):
        CD_3_seq=sequnce[CD_3_seq]
        CD_3_seq=float(len(''.join(CD_3_seq)))
    else:
        CD_3_seq=np.nan
    CD_4D_seq_1=codon_df_4d[(codon_df_4d['CHROM']==window.CHROM) & (codon_df_4d['POS']>=window.BIN_START) & (codon_df_4d['POS']<=window.BIN_END) ]
    CD_4D_seq=list(CD_4D_seq_1['POS']-1)
    if len(CD_4D_seq) != 0 and max(CD_4D_seq) < len(sequnce):
        CD_4D_seq=sequnce[CD_4D_seq]
        CD_4D_seq=float(len(''.join(CD_4D_seq)))
    else:
        CD_4D_seq=np.nan
    list_pos=range(int(window.BIN_START-1), int(window.BIN_END-1))
    merge=pandas.concat([introns_seq_1,CD_1_seq_1,CD_2_seq_1,CD_3_seq_1,CD_4D_seq_1],ignore_index=True)
    intergenic=main_list = np.setdiff1d(list_pos,list(merge.POS))
    intergenic_df=pandas.DataFrame()
    intergenic_df['CHROM']=[window.CHROM]*len(intergenic)
    intergenic_df['POS']=list(intergenic)
    intergenic_final=pandas.concat([intergenic_final, intergenic_df], ignore_index=True)
    intergenic_seq=list(intergenic-1)
    if len(intergenic_seq) != 0 and max(intergenic_seq) < len(sequnce):
        intergenic_seq=sequnce[intergenic_seq]
        intergenic_seq=float(len(''.join(intergenic_seq)))
    else:
        intergenic_seq=np.nan
    #print str(gc_fraction)+ '   '+ window.CHROM + '     '+ str(window.BIN_START) + '        '+ str(window.BIN_END)
    table_windows_2['reference'].loc[window_index]=fraction
    table_windows_2['intron'].loc[window_index]=introns_seq
    table_windows_2['codon_1'].loc[window_index]=CD_1_seq
    table_windows_2['codon_2'].loc[window_index]=CD_2_seq
    table_windows_2['codon_3'].loc[window_index]=CD_3_seq
    table_windows_2['codon_4D'].loc[window_index]=CD_4D_seq
    table_windows_2['intergenic'].loc[window_index]=intergenic_seq
    print window_index
    sys.stdout.flush()



def callcualte_windowed_pi_from_site_pi(population_pi, dataframe, sites_list, output):
	win_all_sites=dataframe[['CHROM','BIN_START','BIN_END',sites_list]]
	win_all_sites=win_all_sites[win_all_sites[sites_list]>0.0]
	win_sites=population_pi
	df_merge=pandas.DataFrame(columns=['CHROM','POS','PI','BIN_START','BIN_END','Sites'])
	l=list(set(win_all_sites.CHROM))
	l.sort()
	for i in l:
		chrom_set=win_all_sites[win_all_sites['CHROM']==i]
		chrom_sites=win_sites[win_sites['CHROM']==i]
		new = pandas.merge(chrom_set, chrom_sites, on=['CHROM'],how='inner')
		new = new.query('BIN_START <= POS and BIN_END >= POS')
		df_merge=pandas.concat([new, df_merge], ignore_index=True)
	out=pandas.DataFrame({output+'_Window_pi' :df_merge.groupby(['CHROM','BIN_START','BIN_END',sites_list], as_index=True)['PI'].sum()}).reset_index()
	out[output+'_Window_pi']=out[output+'_Window_pi']/out[sites_list]
	return out



West =sys.argv[2]

West_site_PI_table =pandas.read_table(West)


introns=introns.merge(West_site_PI_table, on=['CHROM','POS'])
codon_df_1=codon_df_1.merge(West_site_PI_table, on=['CHROM','POS'])
codon_df_2=codon_df_2.merge(West_site_PI_table, on=['CHROM','POS'])
codon_df_3=codon_df_3.merge(West_site_PI_table, on=['CHROM','POS'])
codon_df_4d=codon_df_4d.merge(West_site_PI_table, on=['CHROM','POS'])
intergenic_final=intergenic_final.merge(West_site_PI_table, on=['CHROM','POS'])



#calculate gene wise pi for east and west
introns=callcualte_windowed_pi_from_site_pi(introns, table_windows_2, 'intron', 'intron_West_')
codon_df_1=callcualte_windowed_pi_from_site_pi(codon_df_1, table_windows_2, 'codon_1', 'codon_1_West_')
codon_df_2=callcualte_windowed_pi_from_site_pi(codon_df_2, table_windows_2, 'codon_2', 'codon_2_West_')
codon_df_3=callcualte_windowed_pi_from_site_pi(codon_df_3, table_windows_2, 'codon_3', 'codon_3_West_')
codon_df_4d=callcualte_windowed_pi_from_site_pi(codon_df_4d, table_windows_2, 'codon_4D', 'codon_4D_West_')
intergenic_final=callcualte_windowed_pi_from_site_pi(intergenic_final, table_windows_2, 'intergenic', 'intergenic_final_West_')


introns.to_csv('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/introns_pi_'+str(sys.argv[3])+'_'+str(sys.argv[1])+'.csv')
codon_df_1.to_csv('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/codon_df_1_pi_'+str(sys.argv[3])+'_'+str(sys.argv[1])+'.csv')
codon_df_2.to_csv('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/codon_df_2_pi_'+str(sys.argv[3])+'_'+str(sys.argv[1])+'.csv')
codon_df_3.to_csv('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/codon_df_3_pi_'+str(sys.argv[3])+'_'+str(sys.argv[1])+'.csv')
codon_df_4d.to_csv('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/codon_df_4d_pi_'+str(sys.argv[3])+'_'+str(sys.argv[1])+'.csv')
intergenic_final.to_csv('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/intergenic_final_pi_'+str(sys.argv[3])+'_'+str(sys.argv[1])+'.csv')



