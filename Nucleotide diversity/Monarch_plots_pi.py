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




introns_East=pandas.DataFrame()
codon_df_1_East=pandas.DataFrame()
codon_df_2_East=pandas.DataFrame()
codon_df_3_East=pandas.DataFrame()
codon_df_4d_East=pandas.DataFrame()
intergenic_final_East=pandas.DataFrame()


for i in [0, 1, 2, 4, 5, 6, 7, 8, 9]:
	introns=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/introns_pi_East_'+str(i)+'.csv')
	codon_df_1=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/codon_df_1_pi_East_'+str(i)+'.csv')
	codon_df_2=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/codon_df_2_pi_East_'+str(i)+'.csv')
	codon_df_3=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/codon_df_3_pi_East_'+str(i)+'.csv')
	codon_df_4d=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/codon_df_4d_pi_East_'+str(i)+'.csv')
	intergenic_final=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/intergenic_final_pi_East_'+str(i)+'.csv')
	introns_East=pandas.concat([introns_East, introns], ignore_index=True)
	codon_df_1_East=pandas.concat([codon_df_1_East, codon_df_1], ignore_index=True)
	codon_df_2_East=pandas.concat([codon_df_2_East, codon_df_2], ignore_index=True)
	codon_df_3_East=pandas.concat([codon_df_3_East, codon_df_3], ignore_index=True)
	codon_df_4d_East=pandas.concat([codon_df_4d_East, codon_df_4d], ignore_index=True)
	intergenic_final_East=pandas.concat([intergenic_final_East, intergenic_final], ignore_index=True)


introns_East=introns_East.rename(columns={'intron_West__Window_pi': 'intron_East_Window_pi'})
codon_df_1_East=codon_df_1_East.rename(columns={'codon_1_West__Window_pi': 'codon_1_East_Window_pi'})
codon_df_2_East=codon_df_2_East.rename(columns={'codon_2_West__Window_pi': 'codon_2_East_Window_pi'})
codon_df_3_East=codon_df_3_East.rename(columns={'codon_3_West__Window_pi': 'codon_3_East_Window_pi'})
codon_df_4d_East=codon_df_4d_East.rename(columns={'codon_4D_West__Window_pi': 'codon_4D_East_Window_pi'})
intergenic_final_East=intergenic_final_East.rename(columns={'intergenic_final_West__Window_pi': 'intergenic_final_East_Window_pi'})



introns_West=pandas.DataFrame()
codon_df_1_West=pandas.DataFrame()
codon_df_2_West=pandas.DataFrame()
codon_df_3_West=pandas.DataFrame()
codon_df_4d_West=pandas.DataFrame()
intergenic_final_West=pandas.DataFrame()


for i in [0, 1, 2, 4, 5, 6, 7, 8, 9]:
	introns=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/introns_pi_West'+str(i)+'.csv')
	codon_df_1=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/codon_df_1_pi_West'+str(i)+'.csv')
	codon_df_2=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/codon_df_2_pi_West'+str(i)+'.csv')
	codon_df_3=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/codon_df_3_pi_West'+str(i)+'.csv')
	codon_df_4d=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/codon_df_4d_pi_West'+str(i)+'.csv')
	intergenic_final=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/intergenic_final_pi_West'+str(i)+'.csv')
	introns_West=pandas.concat([introns_West, introns], ignore_index=True)
	codon_df_1_West=pandas.concat([codon_df_1_West, codon_df_1], ignore_index=True)
	codon_df_2_West=pandas.concat([codon_df_2_West, codon_df_2], ignore_index=True)
	codon_df_3_West=pandas.concat([codon_df_3_West, codon_df_3], ignore_index=True)
	codon_df_4d_West=pandas.concat([codon_df_4d_West, codon_df_4d], ignore_index=True)
	intergenic_final_West=pandas.concat([intergenic_final_West, intergenic_final], ignore_index=True)


introns_West=introns_West.rename(columns={'intron_West__Window_pi': 'intron_West_Window_pi'})
codon_df_1_West=codon_df_1_West.rename(columns={'codon_1_West__Window_pi': 'codon_1_West_Window_pi'})
codon_df_2_West=codon_df_2_West.rename(columns={'codon_2_West__Window_pi': 'codon_2_West_Window_pi'})
codon_df_3_West=codon_df_3_West.rename(columns={'codon_3_West__Window_pi': 'codon_3_West_Window_pi'})
codon_df_4d_West=codon_df_4d_West.rename(columns={'codon_4D_West__Window_pi': 'codon_4D_West_Window_pi'})
intergenic_final_West=intergenic_final_West.rename(columns={'intergenic_final_West__Window_pi': 'intergenic_final_West_Window_pi'})



all_info=introns_East.merge(codon_df_1_East, on=['CHROM','BIN_START','BIN_END'], how='outer')
all_info=all_info.merge(codon_df_2_East, on=['CHROM','BIN_START','BIN_END'], how='outer')
all_info=all_info.merge(codon_df_3_East, on=['CHROM','BIN_START','BIN_END'], how='outer')
all_info=all_info.merge(codon_df_4d_East, on=['CHROM','BIN_START','BIN_END'], how='outer')
all_info=all_info.merge(intergenic_final_East, on=['CHROM','BIN_START','BIN_END'], how='outer')
#all_info=all_info[['CHROM','BIN_START','BIN_END','intron','intron_East_Window_pi','codon_1','codon_1_East_Window_pi','codon_2','codon_2_East_Window_pi','codon_3','codon_3_East_Window_pi','codon_4D','codon_4D_East_Window_pi','intergenic','intergenic_final_East_Window_pi']]
all_info=all_info.merge(introns_West, on=['CHROM','BIN_START','BIN_END'], how='outer')
all_info=all_info.merge(codon_df_1_West, on=['CHROM','BIN_START','BIN_END'], how='outer')
all_info=all_info.merge(codon_df_2_West, on=['CHROM','BIN_START','BIN_END'], how='outer')
all_info=all_info.merge(codon_df_3_West, on=['CHROM','BIN_START','BIN_END'], how='outer')
all_info=all_info.merge(codon_df_4d_West, on=['CHROM','BIN_START','BIN_END'], how='outer')
all_info=all_info.merge(intergenic_final_West, on=['CHROM','BIN_START','BIN_END'], how='outer')

all_info=all_info[['CHROM','BIN_START','BIN_END','intron_East_Window_pi','codon_1_East_Window_pi','codon_2_East_Window_pi','codon_3_East_Window_pi','codon_4D_East_Window_pi','intergenic_final_East_Window_pi','intron_West_Window_pi','codon_1_West_Window_pi','codon_2_West_Window_pi','codon_3_West_Window_pi','codon_4D_West_Window_pi','intergenic_final_West_Window_pi']]

Fst_Dxy_windows=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/Final_db.csv')

Fst_Dxy_windows=Fst_Dxy_windows[(Fst_Dxy_windows.East_N_sites > 3000) & (Fst_Dxy_windows.West_N_sites > 3000)]

all_info=Fst_Dxy_windows.merge(all_info, on=['CHROM','BIN_START','BIN_END'])
all_info=all_info[['CHROM','BIN_START','BIN_END','intron_East_Window_pi','codon_1_East_Window_pi','codon_2_East_Window_pi','codon_3_East_Window_pi','codon_4D_East_Window_pi','intergenic_final_East_Window_pi','intron_West_Window_pi','codon_1_West_Window_pi','codon_2_West_Window_pi','codon_3_West_Window_pi','codon_4D_West_Window_pi','intergenic_final_West_Window_pi','East_final_pi','West_final_pi']]


plotstuff=all_info[['intron_East_Window_pi','codon_1_East_Window_pi','codon_2_East_Window_pi','codon_3_East_Window_pi','codon_4D_East_Window_pi','intergenic_final_East_Window_pi','intron_West_Window_pi','codon_1_West_Window_pi','codon_2_West_Window_pi','codon_3_West_Window_pi','codon_4D_West_Window_pi','intergenic_final_West_Window_pi','East_final_pi','West_final_pi']]


violin_df=pandas.DataFrame()
violin_df['PI']=list(plotstuff['intron_East_Window_pi'])+list(plotstuff['codon_1_East_Window_pi'])+list(plotstuff['codon_2_East_Window_pi'])+list(plotstuff['codon_3_East_Window_pi'])+list(plotstuff['codon_4D_East_Window_pi'])+list(plotstuff['intergenic_final_East_Window_pi'])+list(plotstuff['East_final_pi'])+list(plotstuff['intron_West_Window_pi'])+list(plotstuff['codon_1_West_Window_pi'])+list(plotstuff['codon_2_West_Window_pi'])+list(plotstuff['codon_3_West_Window_pi'])+list(plotstuff['codon_4D_West_Window_pi'])+list(plotstuff['intergenic_final_West_Window_pi'])+list(plotstuff['West_final_pi'])
violin_df['CAT']=list(['intron']*len(list(plotstuff['intron_East_Window_pi'])))+list(['codon_1']*len(list(plotstuff['codon_1_East_Window_pi'])))+list(['codon_2']*len(list(plotstuff['codon_2_East_Window_pi'])))+list(['codon_3']*len(list(plotstuff['codon_3_East_Window_pi'])))+list(['codon_4D']*len(list(plotstuff['codon_4D_East_Window_pi'])))+list(['intergenic_final']*len(list(plotstuff['intergenic_final_East_Window_pi'])))+list(['all']*len(list(plotstuff['East_final_pi'])))+list(['intron']*len(list(plotstuff['intron_West_Window_pi'])))+list(['codon_1']*len(list(plotstuff['codon_1_West_Window_pi'])))+list(['codon_2']*len(list(plotstuff['codon_2_West_Window_pi'])))+list(['codon_3']*len(list(plotstuff['codon_3_West_Window_pi'])))+list(['codon_4D']*len(list(plotstuff['codon_4D_West_Window_pi'])))+list(['intergenic_final']*len(list(plotstuff['intergenic_final_West_Window_pi'])))+list(['all']*len(list(plotstuff['West_final_pi'])))
violin_df['POP']=list(['East']*len(list(plotstuff['intron_East_Window_pi'])))+list(['East']*len(list(plotstuff['codon_1_East_Window_pi'])))+list(['East']*len(list(plotstuff['codon_2_East_Window_pi'])))+list(['East']*len(list(plotstuff['codon_3_East_Window_pi'])))+list(['East']*len(list(plotstuff['codon_4D_East_Window_pi'])))+list(['East']*len(list(plotstuff['intergenic_final_East_Window_pi'])))+list(['East']*len(list(plotstuff['East_final_pi'])))+list(['West']*len(list(plotstuff['intron_West_Window_pi'])))+list(['West']*len(list(plotstuff['codon_1_West_Window_pi'])))+list(['West']*len(list(plotstuff['codon_2_West_Window_pi'])))+list(['West']*len(list(plotstuff['codon_3_West_Window_pi'])))+list(['West']*len(list(plotstuff['codon_4D_West_Window_pi'])))+list(['West']*len(list(plotstuff['intergenic_final_West_Window_pi'])))+list(['West']*len(list(plotstuff['West_final_pi'])))

ax = sns.boxplot(x="CAT", y="PI", hue="POP",data=violin_df,order=['all','intergenic_final','intron','codon_1','codon_2','codon_3','codon_4D'], palette=sns.color_palette(['#737070','#d90d0d']),showfliers=False)




