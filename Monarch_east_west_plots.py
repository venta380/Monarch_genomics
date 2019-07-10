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
import re

sns.set(font_scale=1.5)
sns.set_style("whitegrid", {'axes.grid' : False})

pwd='/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp'
os.chdir(pwd)
def join_raw_data_base(lists):
        for i in range(1,len(lists)):
                j=i-1
                if i == 1:
                        new_2=pandas.merge(lists[j], lists[i], on=['CHROM','BIN_START'], how='inner')
                else:
                        new_2=pandas.merge(nextone, lists[i], on=['CHROM','BIN_START'], how='inner')
                nextone=new_2
        return nextone


def corr_pvalue(df):
        from scipy.stats import pearsonr
        import numpy as np
        import pandas as pd
        numeric_df = df.dropna()._get_numeric_data()
        cols = numeric_df.columns
        mat = numeric_df.values
        arr = np.zeros((len(cols),len(cols)), dtype=object)
        for xi, x in enumerate(mat.T):
                for yi, y in enumerate(mat.T[xi:]):
                    arr[xi, yi+xi] = map(lambda _: round(_,3), pearsonr(x,y))
                    arr[yi+xi, xi] = arr[xi, yi+xi]
        
        return pandas.DataFrame(arr, index=cols, columns=cols)

def corr_plot(df):
        from pandas.plotting import scatter_matrix
        axes= scatter_matrix(df,diagonal='none')
        corr = df.corr(method='pearson').as_matrix()
        mask = np.zeros_like(corr, dtype=np.bool)
        mask[np.tril_indices_from(mask)] = True
        #for i, j in zip(*plt.np.tril_indices_from(axes, k=1)):
        #    axes[i, j].annotate("%.3f" %corr[i,j], (0.1, 0.8), xycoords='axes fraction', ha='center', va='center')
        
        for i, j in zip(*plt.np.triu_indices_from(axes, k=1)):
            axes[i, j].set_visible(False)
        
        for i, j in zip(*plt.np.diag_indices_from(axes)):
            axes[i, j].set_visible(False)

        for i, j in zip(*plt.np.tril_indices_from(axes, k=1)):
            axes[i, j].set(xlim=(0,0.005))
            axes[i, j].set(ylim=(0,0.005))
        plt.show()



def genes_in_outliers(data):
        head_stuff=['scaffold','source','feature','start','end','.','orentation','..','infor']
        annotation=pandas.read_table("/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_genes.gff",skiprows=0, names=head_stuff)
        annotation=annotation[['scaffold','source','feature','start','end','infor','orentation']]
        annotation=annotation[annotation.feature=='gene']
        list_genes=[]
        list_CHROM=[]
        list_BIN_START=[]
        list_BIN_END=[]
        for window_index, window_1 in data.iterrows():
                window=annotation[(annotation['scaffold'] == window_1.CHROM)]
                window=window[((window['start'].isin(range(int(window_1.BIN_START), int(window_1.BIN_END))))) | ((window['end'].isin(range(int(window_1.BIN_START), int(window_1.BIN_END)))))]
                #window=window.append(window[(window['start']<=int(window_1.BIN_START)) & (window['end']>=int(window_1.BIN_END))])
                #print window.empty
                if window.empty != 'False':
                        #print window
                        #print window_1.BIN_START
                        #print window_1.BIN_END
                        for i_index, i in window.iterrows():
                                list_genes.append(str(i['infor']))
                                list_CHROM.append(str(window_1.CHROM))
                                list_BIN_START.append(str(window_1.BIN_START))
                                list_BIN_END.append(str(window_1.BIN_END))
        list_genes_=[i.split(';')[1][5:] for i in list_genes]
        #f=open('data_Z_FST_16.fasta', 'w+')
        file_name='ZFST_'+list(data)[2].split('_')[1]+'_genes_in_outliers.csv'
        d = {'gene': list_genes_ , 'CHROM': list_CHROM , 'BIN_START': list_BIN_START , 'BIN_END': list_BIN_END}
        df = pandas.DataFrame(data=d)
        #df.to_csv(file_name, header=False, index=False,sep = " ")
        return df 


windows=pandas.read_csv("windows.csv")
windows_N_sites=pandas.read_csv('./bam_covrage_files/East_West_N_sites')

FST_windows=pandas.read_table("Fst_window_10000.windowed.weir.fst")
FST_windows=pandas.merge(windows_N_sites, FST_windows, on=['CHROM','BIN_START','BIN_END'], how='inner')

FST_windows=FST_windows[FST_windows.East_West_N_sites >= 1000]

FST_windows['Scaled_FST']=(FST_windows['MEAN_FST']*10000)/FST_windows['East_West_N_sites']
FST_windows['Scaled_FST']=FST_windows['Scaled_FST'].clip(lower=0)


Td_windows_East=pandas.read_table("East.Td.sites.Tajima.D")
Td_windows_East=Td_windows_East.rename(index=str, columns={'TajimaD': 'East_TajimaD' })
Td_windows_West=pandas.read_table("random_west.Td.sites.Tajima.D")
Td_windows_West=Td_windows_West.rename(index=str, columns={'TajimaD': 'West_TajimaD' })
Td_windows_East['BIN_START']=Td_windows_East.BIN_START+1
Td_windows_West['BIN_START']=Td_windows_West.BIN_START+1
Td_windows_East['BIN_END']=Td_windows_East.BIN_START+9999
Td_windows_West['BIN_END']=Td_windows_West.BIN_START+9999

East_allele_frquencey ="/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/East.frq"
West_allele_frquencey ="/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/random_west.frq"

FST_windows[FST_windows['N_VARIANTS']>100]['MEAN_FST'].mean()
FST_windows=pandas.merge(FST_windows, Td_windows_East, on=['CHROM', 'BIN_START', 'BIN_END'], how='inner')
FST_windows=pandas.merge(FST_windows, Td_windows_West, on=['CHROM', 'BIN_START', 'BIN_END'], how='inner')


def callcualte_windowed_FST_from_site_FST(population_pi, dataframe, sites_list, output):
	win_all_sites=dataframe[['CHROM','BIN_START','BIN_END',sites_list]]
	win_sites=population_pi
	win_sites['BIN_START']=(np.floor(win_sites['POS']/10000)*10000)+1
	win_sites['BIN_END']=win_sites['BIN_START']+(10000-1)
	df_merge = pandas.merge(win_sites, win_all_sites, on=['CHROM','BIN_START','BIN_END'],how='inner')
	#df_merge = df_merge.query('BIN_START <= POS and BIN_END >= POS')
	out=pandas.DataFrame({output :df_merge.groupby(['CHROM','BIN_START','BIN_END',sites_list], as_index=True)['WEIR_AND_COCKERHAM_FST'].sum()}).reset_index()
	out[output]=out[output]/out[sites_list]
	return out


######load sites covered in population comparisions
Big_Sur_VS_Carpinteria_N_sites=pandas.read_csv('./bam_covrage_files/Big_Sur_VS_Carpinteria_N_sites')
Big_Sur_VS_Oceano_N_sites=pandas.read_csv('./bam_covrage_files/Big_Sur_VS_Oceano_N_sites')
Carpinteria_VS_Oceano_N_sites=pandas.read_csv('./bam_covrage_files/Carpinteria_VS_Oceano_N_sites')
Carpinteria_VS_East_N_sites=pandas.read_csv('./bam_covrage_files/Carpinteria_VS_East_N_sites')
Oceano_VS_East_N_sites=pandas.read_csv('./bam_covrage_files/Oceano_VS_East_N_sites')
big_sur_VS_East_N_sites=pandas.read_csv('./bam_covrage_files/big_sur_VS_East_N_sites')

FST_windows=pandas.merge(FST_windows, Big_Sur_VS_Carpinteria_N_sites, on=['CHROM', 'BIN_START', 'BIN_END'], how='inner')
FST_windows=pandas.merge(FST_windows, Big_Sur_VS_Oceano_N_sites, on=['CHROM', 'BIN_START', 'BIN_END'], how='inner')
FST_windows=pandas.merge(FST_windows, Carpinteria_VS_Oceano_N_sites, on=['CHROM', 'BIN_START', 'BIN_END'], how='inner')
FST_windows=pandas.merge(FST_windows, Carpinteria_VS_East_N_sites, on=['CHROM', 'BIN_START', 'BIN_END'], how='inner')
FST_windows=pandas.merge(FST_windows, Oceano_VS_East_N_sites, on=['CHROM', 'BIN_START', 'BIN_END'], how='inner')
FST_windows=pandas.merge(FST_windows, big_sur_VS_East_N_sites, on=['CHROM', 'BIN_START', 'BIN_END'], how='inner')


Fst_sites=pandas.read_table('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/Fst_window_site.weir.fst')
Fst_sites['WEIR_AND_COCKERHAM_FST']=Fst_sites['WEIR_AND_COCKERHAM_FST'].clip(lower=0)

new=callcualte_windowed_FST_from_site_FST(Fst_sites, FST_windows, 'East_West_N_sites', 'WEIR_AND_COCKERHAM_MEAN_FST')

FST_windows=pandas.merge(FST_windows, new, on=['CHROM', 'BIN_START', 'BIN_END'], how='inner')

#######################
Fst_sites_Carpinteria_VS_Oceano=pandas.read_table('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/WEST_pop_Carpinteria_VS_Oceano.weir.fst')
Fst_sites_Carpinteria_VS_Oceano['WEIR_AND_COCKERHAM_FST']=Fst_sites_Carpinteria_VS_Oceano['WEIR_AND_COCKERHAM_FST'].clip(lower=0)

Carpinteria_VS_Oceano=callcualte_windowed_FST_from_site_FST(Fst_sites_Carpinteria_VS_Oceano, FST_windows, 'Carpinteria_VS_Oceano_N_sites', 'Carpinteria_VS_Oceano_MEAN_FST')

Carpinteria_VS_Oceano=pandas.merge(FST_windows[['CHROM', 'BIN_START', 'BIN_END']], Carpinteria_VS_Oceano, on=['CHROM', 'BIN_START', 'BIN_END'], how='inner')
#0.000908 ± 0.000381
#######################

#######################
Fst_sites_Big_Sur_VS_Oceano=pandas.read_table('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/WEST_pop_Big_Sur_VS_Oceano.weir.fst')
Fst_sites_Big_Sur_VS_Oceano['WEIR_AND_COCKERHAM_FST']=Fst_sites_Big_Sur_VS_Oceano['WEIR_AND_COCKERHAM_FST'].clip(lower=0)

Big_Sur_VS_Oceano=callcualte_windowed_FST_from_site_FST(Fst_sites_Big_Sur_VS_Oceano, FST_windows, 'Big_Sur_VS_Oceano_N_sites', 'Big_Sur_VS_Oceano_MEAN_FST')

Big_Sur_VS_Oceano=pandas.merge(FST_windows[['CHROM', 'BIN_START', 'BIN_END']], Big_Sur_VS_Oceano, on=['CHROM', 'BIN_START', 'BIN_END'], how='inner')
#0.001237 ± 0.000516
#######################

#######################
Fst_sites_Big_Sur_VS_Carpinteria=pandas.read_table('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/WEST_pop_Big_Sur_VS_Carpinteria.weir.fst')
Fst_sites_Big_Sur_VS_Carpinteria['WEIR_AND_COCKERHAM_FST']=Fst_sites_Big_Sur_VS_Carpinteria['WEIR_AND_COCKERHAM_FST'].clip(lower=0)

Big_Sur_VS_Carpinteria=callcualte_windowed_FST_from_site_FST(Fst_sites_Big_Sur_VS_Carpinteria, FST_windows, 'Big_Sur_VS_Carpinteria_N_sites', 'Big_Sur_VS_Carpinteria_MEAN_FST')

Big_Sur_VS_Carpinteria=pandas.merge(FST_windows[['CHROM', 'BIN_START', 'BIN_END']], Big_Sur_VS_Carpinteria, on=['CHROM', 'BIN_START', 'BIN_END'], how='inner')
#0.001332 ± 0.000544
#######################
#######################
Fst_sites_Carpinteria_VS_East=pandas.read_table('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/WEST_pop_Carpinteria_VS_East.weir.fst')
Fst_sites_Carpinteria_VS_East['WEIR_AND_COCKERHAM_FST']=Fst_sites_Carpinteria_VS_East['WEIR_AND_COCKERHAM_FST'].clip(lower=0)

Carpinteria_VS_East=callcualte_windowed_FST_from_site_FST(Fst_sites_Carpinteria_VS_East, FST_windows, 'Carpinteria_VS_East_N_sites', 'Carpinteria_VS_East_MEAN_FST')

Carpinteria_VS_East=pandas.merge(FST_windows[['CHROM', 'BIN_START', 'BIN_END']], Carpinteria_VS_East, on=['CHROM', 'BIN_START', 'BIN_END'], how='inner')
#######################
#######################
Fst_sites_Oceano_VS_East=pandas.read_table('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/WEST_pop_Oceano_VS_East.weir.fst')
Fst_sites_Oceano_VS_East['WEIR_AND_COCKERHAM_FST']=Fst_sites_Oceano_VS_East['WEIR_AND_COCKERHAM_FST'].clip(lower=0)

Oceano_VS_East=callcualte_windowed_FST_from_site_FST(Fst_sites_Oceano_VS_East, FST_windows, 'Oceano_VS_East_N_sites', 'Oceano_VS_East_MEAN_FST')

Oceano_VS_East=pandas.merge(FST_windows[['CHROM', 'BIN_START', 'BIN_END']], Oceano_VS_East, on=['CHROM', 'BIN_START', 'BIN_END'], how='inner')
#######################
#######################
Fst_sites_big_sur_VS_East=pandas.read_table('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/WEST_pop_big_sur_VS_East.weir.fst')
Fst_sites_big_sur_VS_East['WEIR_AND_COCKERHAM_FST']=Fst_sites_big_sur_VS_East['WEIR_AND_COCKERHAM_FST'].clip(lower=0)

big_sur_VS_East=callcualte_windowed_FST_from_site_FST(Fst_sites_big_sur_VS_East, FST_windows, 'big_sur_VS_East_N_sites', 'big_sur_VS_East_MEAN_FST')

big_sur_VS_East=pandas.merge(FST_windows[['CHROM', 'BIN_START', 'BIN_END']], big_sur_VS_East, on=['CHROM', 'BIN_START', 'BIN_END'], how='inner')
#######################

population_comparisions=pandas.merge(Carpinteria_VS_Oceano, Big_Sur_VS_Oceano,on=['CHROM', 'BIN_START', 'BIN_END'], how='inner')
population_comparisions=pandas.merge(population_comparisions, Big_Sur_VS_Carpinteria,on=['CHROM', 'BIN_START', 'BIN_END'], how='inner')
population_comparisions=pandas.merge(population_comparisions, Carpinteria_VS_East,on=['CHROM', 'BIN_START', 'BIN_END'], how='inner')
population_comparisions=pandas.merge(population_comparisions, Oceano_VS_East,on=['CHROM', 'BIN_START', 'BIN_END'], how='inner')
population_comparisions=pandas.merge(population_comparisions, big_sur_VS_East,on=['CHROM', 'BIN_START', 'BIN_END'], how='inner')

#sns.boxplot(data=population_comparisions[['Carpinteria_VS_Oceano_MEAN_FST','Big_Sur_VS_Oceano_MEAN_FST','Big_Sur_VS_Carpinteria_MEAN_FST','Carpinteria_VS_East_MEAN_FST','Oceano_VS_East_MEAN_FST','big_sur_VS_East_MEAN_FST']],showfliers=False, palette=sns.color_palette(['#CA9261','#F8B1E2','#949494','#EEE030','#56B4E8','#0173B0']))


West_pop=pandas.merge(FST_windows[['CHROM', 'BIN_START', 'BIN_END', 'WEIR_AND_COCKERHAM_MEAN_FST']], Carpinteria_VS_Oceano, on=['CHROM', 'BIN_START', 'BIN_END'], how='inner')
West_pop=pandas.merge(West_pop, Big_Sur_VS_Oceano, on=['CHROM', 'BIN_START', 'BIN_END'], how='inner')
West_pop=pandas.merge(West_pop, Big_Sur_VS_Carpinteria, on=['CHROM', 'BIN_START', 'BIN_END'], how='inner')

West_pop.mean()
West_pop=West_pop[['CHROM', 'BIN_START', 'BIN_END', 'Carpinteria_VS_Oceano_MEAN_FST', 'Big_Sur_VS_Oceano_MEAN_FST', 'Big_Sur_VS_Carpinteria_MEAN_FST']]
#WEIR_AND_COCKERHAM_MEAN_FST             0.001041
#Carpinteria_VS_Oceano_MEAN_FST          0.000908
#Big_Sur_VS_Oceano_MEAN_FST              0.001237
#Big_Sur_VS_Carpinteria_MEAN_FST         0.001332


Chrom=pandas.read_csv("/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/Chromosomes_final_2.txt",names=['CHROM','Auto_Allo','Chr_number'], sep=',',header=None)



Dxy_windows=pandas.read_csv("East_West_Dxy_DB")
Pi_windows=pandas.read_csv("final_DB_pi.csv")
Pi_windows=Pi_windows[['CHROM', 'BIN_START', 'BIN_END', 'East_final_pi', 'West_final_pi']]


Fst_Dxy_windows=pandas.merge(FST_windows, Dxy_windows, on=['CHROM', 'BIN_START', 'BIN_END'], how='inner')
Fst_Dxy_windows=pandas.merge(Fst_Dxy_windows, West_pop, on=['CHROM', 'BIN_START', 'BIN_END'], how='inner')

#Fst_Dxy_windows=Fst_Dxy_windows.rename(columns={'East_West_N_sites_x': 'East_West_N_sites'})
Fst_Dxy_windows=pandas.merge(Fst_Dxy_windows, Pi_windows, on=['CHROM', 'BIN_START', 'BIN_END'], how='inner')


Fst_Dxy_windows.loc[(Fst_Dxy_windows['CHROM']=='DPSCF300001')&(Fst_Dxy_windows['BIN_START']<5820000), 'CHROM']='DPSCF300001-1'
Fst_Dxy_windows.loc[(Fst_Dxy_windows['CHROM']=='DPSCF300001')&(Fst_Dxy_windows['BIN_START']>5820000), 'CHROM']='DPSCF300001-2'

Fst_Dxy_windows.loc[(Fst_Dxy_windows['CHROM']=='DPSCF300028')&(Fst_Dxy_windows['BIN_START']<760000), 'CHROM']='DPSCF300028-1'
Fst_Dxy_windows.loc[(Fst_Dxy_windows['CHROM']=='DPSCF300028')&(Fst_Dxy_windows['BIN_START']>1805000), 'CHROM']='DPSCF300028-3'
Fst_Dxy_windows.loc[(Fst_Dxy_windows['CHROM']=='DPSCF300028')&(Fst_Dxy_windows['BIN_START']>760000)&(Fst_Dxy_windows['BIN_START']<1805000), 'CHROM']='DPSCF300028-2'

Fst_Dxy_windows.loc[(Fst_Dxy_windows['CHROM']=='DPSCF300044')&(Fst_Dxy_windows['BIN_START']<290000), 'CHROM']='DPSCF300044-1'
Fst_Dxy_windows.loc[(Fst_Dxy_windows['CHROM']=='DPSCF300044')&(Fst_Dxy_windows['BIN_START']>290000), 'CHROM']='DPSCF300044-2'

Fst_Dxy_windows=pandas.merge(Fst_Dxy_windows, Chrom, on=['CHROM'], how='inner')

Fst_Dxy_windows.loc[(Fst_Dxy_windows['CHROM']=='DPSCF300001-1')&(Fst_Dxy_windows['BIN_START']>3827349)&(Fst_Dxy_windows['BIN_END']<3882816), 'Chr_number']='Z_neo'
Fst_Dxy_windows.loc[(Fst_Dxy_windows['CHROM']=='DPSCF300001-1')&(Fst_Dxy_windows['BIN_START']>3827349)&(Fst_Dxy_windows['BIN_END']<3882816), 'Auto_Allo']='Z_neo'
Fst_Dxy_windows.loc[(Fst_Dxy_windows['CHROM']=='DPSCF300001-1')&(Fst_Dxy_windows['BIN_START']>3827349)&(Fst_Dxy_windows['BIN_END']<3882816), 'CHROM']='DPSCF300001-Z_neo'
Fst_Dxy_windows.loc[(Fst_Dxy_windows['Auto_Allo']=='Z_neo'), 'Chr_number']='Z_neo'

Fst_Dxy_windows=Fst_Dxy_windows[['CHROM','BIN_START','BIN_END','East_West_N_sites','West_N_sites','East_N_sites','N_VARIANTS','WEIGHTED_FST','MEAN_FST','Scaled_FST','East_TajimaD','West_TajimaD','Auto_Allo','Chr_number','East_West_dxy','East_West_fst','East_West_fixed','East_West_private_a','East_West_private_b','East_West_shared','East_West_same_sites','East_West_Pi_a','East_West_Pi_b','East_final_pi','West_final_pi','WEIR_AND_COCKERHAM_MEAN_FST', 'Carpinteria_VS_Oceano_MEAN_FST', 'Big_Sur_VS_Oceano_MEAN_FST', 'Big_Sur_VS_Carpinteria_MEAN_FST']]



CHROM=["1",'Z_neo',"2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","Un_Assigned"]
CHROM_dict={}

for i in CHROM:
	CH=list(set(Fst_Dxy_windows[Fst_Dxy_windows['Chr_number']==i]['CHROM']))
	CH.sort()
	CHROM_dict[i]=CH


def manhattan_plot_fst_dxy_pi(sites_filter_0, fst, dxy,pi,taj,fixed,chromosomes_dict,colour):
	species_colour=sns.color_palette("Set2", 8)[:3]
	chromosome=["1","Z_neo","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","Un_Assigned"]
	pi_selection=pi
	fst_selection=fst+['CHROM']
	dxy_selection=dxy+['CHROM']
	taj_selection=taj
	fst_set=[]
	dxy_set=[]
	fixed_set=[]
	pi_set=[]
	taj_set=[]
	chrom={}
	for j in chromosome:
		chrom[j]=0
		for i in chromosomes_dict[str(j)]:
			pi_list=sites_filter_0[(sites_filter_0['CHROM']==i)][pi_selection]
			fixed_1=sites_filter_0[(sites_filter_0['CHROM']==i)][fixed]
			chrom[j]+=len(sites_filter_0[(sites_filter_0['CHROM']==i)])
			fst_list=sites_filter_0[(sites_filter_0['CHROM']==i)][fst_selection]
			dxy_list=sites_filter_0[(sites_filter_0['CHROM']==i)][dxy_selection]
			taj_list=sites_filter_0[(sites_filter_0['CHROM']==i)][taj_selection]
			if len(fst_list)>=5:
				fst_=fst_list[fst].rolling(window=5).mean()
				fst_set.append(fst_)
			if len(dxy_list)>=5:
				dxy_=dxy_list[dxy].rolling(window=5).mean()
				dxy_set.append(dxy_)
			if len(pi_list)>=5:
				pi_=pi_list.rolling(window=5).mean()
				pi_set.append(pi_)
				#pi_.mean().plot(fontsize=10,legend=True)
			if len(fixed_1)>=1:
				fixed_=fixed_1.rolling(window=5).mean()
				fixed_set.append(fixed_)
			if len(taj_list)>=5:
				taj_=taj_list.rolling(window=5).mean()
				taj_set.append(taj_)
				#pi_.mean().plot(fontsize=10,legend=True)
	fst_set2=pandas.DataFrame(columns=fst_selection)
	for j in fst_set:
		list_j=[fst_set2,j]
		fst_set2=pandas.concat(list_j,ignore_index=True)
	
		
	dxy_set2=pandas.DataFrame(columns=dxy_selection)
	for n in dxy_set:
		list_n=[dxy_set2,n]
		dxy_set2=pandas.concat(list_n,ignore_index=True)

	fixed_set2=pandas.DataFrame(columns=fixed)
	for f in fixed_set:
		list_f=[fixed_set2,f]
		fixed_set2=pandas.concat(list_f,ignore_index=True)

	pi_set2=pandas.DataFrame(columns=pi_selection)
	for p in pi_set:
		list_p=[pi_set2,p]
		pi_set2=pandas.concat(list_p,ignore_index=True)

	taj_set2=pandas.DataFrame(columns=taj_selection)
	for p in taj_set:
		list_p=[taj_set2,p]
		taj_set2=pandas.concat(list_p,ignore_index=True)

	fig, axes = plt.subplots(nrows=5, sharex=True)
	fig.subplots_adjust(hspace=0.1)
	fixed_set2.plot.area(ax=axes[0],legend=False,color=sns.color_palette(['#1E78B4','#C54E52', '#56A867', '#DF8352']),stacked=False, alpha=0.5).set_ylabel('fixed', fontsize=15)
	fst_set2.plot(ax=axes[1],legend=False,color='#f4b400').set_ylabel('Fst', fontsize=15)
	dxy_set2.plot(ax=axes[2], legend=False,color='#f4b400').set_ylabel('Dxy', fontsize=15)
	pi_set2.plot(ax=axes[3], legend=False ,color=sns.color_palette(['#030303','#d90d0d']), alpha=1.0).set_ylabel(r'$\pi$', fontsize=15)
	taj_set2.plot(ax=axes[4], legend=False ,color=sns.color_palette(['#030303','#d90d0d']), alpha=1.0).set_ylabel('TajimaD', fontsize=15)
	li=[]
	for c in chromosome:
		if c == '1':
			start=0
			end=int(chrom[str(c)])
		else:
			start=end
			end=start+int(chrom[str(c)])
		li.append([start, end])
	for l in li[0::2]:
		axes[0].axvspan(l[0], l[1], facecolor='#dbd9d9', alpha=0.5)
		axes[1].axvspan(l[0], l[1], facecolor='#dbd9d9', alpha=0.5)
		axes[2].axvspan(l[0], l[1], facecolor='#dbd9d9', alpha=0.5)
		axes[3].axvspan(l[0], l[1], facecolor='#dbd9d9', alpha=0.5)
		axes[4].axvspan(l[0], l[1], facecolor='#dbd9d9', alpha=0.5)
	ticks = axes[4].get_xticks()/100
	axes[4].set_xticklabels(ticks)
	axes[4].set_xlabel("Mega bases")
	return fig

######

z=Fst_Dxy_windows[(Fst_Dxy_windows['Auto_Allo'] == 'Z_anc')]
N_z=Fst_Dxy_windows[(Fst_Dxy_windows['Auto_Allo'] == 'Z_neo')]
a=Fst_Dxy_windows[(Fst_Dxy_windows['Auto_Allo'] == 'A')]




###################################fixed_private_shared_same##########################################

#labels = 'Fixed', 'Private A', 'Private B', 'shared'
#list_sites=["East_West_"]
#
#
#
#for i in list_sites:
#	j=[i+'fixed', i+'private_a', i+'private_b', i+'shared']
#	n=list(Fst_Dxy_windows[j].sum())
#	print i
#	print list(Fst_Dxy_windows[j].sum())
#	print sum(list(Fst_Dxy_windows[j].sum()))
#	sizes=[(j/np.nansum(n)) *100   for j in n]
#	explode = (0, 0, 0, 0)
#	fig1, ax1 = plt.subplots()
#	ax1.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%', shadow=True, startangle=90)
#	ax1.axis('equal')
#	ax1.set_title(i[:-1], bbox={'facecolor':'0.8', 'pad':5})
#
#
#plt.show()
#

################Zfst

z=Fst_Dxy_windows[(Fst_Dxy_windows['Auto_Allo'] == 'Z_anc')]
N_z=Fst_Dxy_windows[(Fst_Dxy_windows['Auto_Allo'] == 'Z_neo')]
a=Fst_Dxy_windows[(Fst_Dxy_windows['Auto_Allo'] == 'A')]

z['Z_FST']=(z['WEIR_AND_COCKERHAM_MEAN_FST']-np.nanmean(z['WEIR_AND_COCKERHAM_MEAN_FST']))/np.nanstd(z['WEIR_AND_COCKERHAM_MEAN_FST'])
N_z['Z_FST']=(N_z['WEIR_AND_COCKERHAM_MEAN_FST']-np.nanmean(N_z['WEIR_AND_COCKERHAM_MEAN_FST']))/np.nanstd(N_z['WEIR_AND_COCKERHAM_MEAN_FST'])
a['Z_FST']=(a['WEIR_AND_COCKERHAM_MEAN_FST']-np.nanmean(a['WEIR_AND_COCKERHAM_MEAN_FST']))/np.nanstd(a['WEIR_AND_COCKERHAM_MEAN_FST'])

all_a_Z_FST =a[(a['Z_FST']  > a.Z_FST.dropna().quantile(0.99))]
all_N_z_Z_FST =N_z[(N_z['Z_FST']  > N_z.Z_FST.dropna().quantile(0.99))]
all_z_Z_FST =z[(z['Z_FST']  > z.Z_FST.dropna().quantile(0.99))]

all_data_Z_FST =pandas.concat([all_a_Z_FST ,all_N_z_Z_FST, all_z_Z_FST], ignore_index=True)
all_data_Z_FST =all_data_Z_FST[all_data_Z_FST['East_West_N_sites'] >= 2000]


a_Z_FST =a[(a['Z_FST']  > a.Z_FST.dropna().quantile(0.99))]
N_z_Z_FST =N_z[(N_z['Z_FST']  > N_z.Z_FST.dropna().quantile(0.99))]
z_Z_FST =z[(z['Z_FST']  > z.Z_FST.dropna().quantile(0.99))]

data_Z_FST_windows_outliers =pandas.concat([a_Z_FST ,N_z_Z_FST, z_Z_FST], ignore_index=True)

data_Z_FST_windows_outliers_genes=genes_in_outliers(data_Z_FST_windows_outliers)
genes_from_paper=[ 'DPOGS200868', 'DPOGS201379', 'DPOGS211203', 'DPOGS202675', 'DPOGS215054', 'DPOGS211604']

data_Z_FST_windows_outliers_genes[data_Z_FST_windows_outliers_genes.gene.isin(genes_from_paper)]
################Zfst_Writes_fst

z=Fst_Dxy_windows[(Fst_Dxy_windows['Auto_Allo'] == 'Z_anc')]
N_z=Fst_Dxy_windows[(Fst_Dxy_windows['Auto_Allo'] == 'Z_neo')]
a=Fst_Dxy_windows[(Fst_Dxy_windows['Auto_Allo'] == 'A')]


z['Z_FST']=(z['East_West_fst']-np.nanmean(z['East_West_fst']))/np.nanstd(z['East_West_fst'])
N_z['Z_FST']=(N_z['East_West_fst']-np.nanmean(N_z['East_West_fst']))/np.nanstd(N_z['East_West_fst'])
a['Z_FST']=(a['East_West_fst']-np.nanmean(a['East_West_fst']))/np.nanstd(a['East_West_fst'])

all_a_Z_FST =a[(a['Z_FST']  > a.Z_FST.dropna().quantile(0.99))]
all_N_z_Z_FST =N_z[(N_z['Z_FST']  > N_z.Z_FST.dropna().quantile(0.99))]
all_z_Z_FST =z[(z['Z_FST']  > z.Z_FST.dropna().quantile(0.99))]

all_data_Z_FST_Write =pandas.concat([all_a_Z_FST ,all_N_z_Z_FST, all_z_Z_FST], ignore_index=True)
all_data_Z_FST_Write =all_data_Z_FST[all_data_Z_FST.East_West_N_sites >= 2000]

all_data_Z_FST_merge=pandas.merge(data_Z_FST_windows_outliers, all_data_Z_FST_Write, on=['CHROM','BIN_START','BIN_END'])[['CHROM','BIN_START','BIN_END']]


##################################
all_data_Z_FST_merge.loc[(all_data_Z_FST_merge['CHROM']=='DPSCF300001-1'), 'CHROM']='DPSCF300001'
all_data_Z_FST_merge.loc[(all_data_Z_FST_merge['CHROM']=='DPSCF300001-2'), 'CHROM']='DPSCF300001'
all_data_Z_FST_merge.loc[(all_data_Z_FST_merge['CHROM']=='DPSCF300028-1'), 'CHROM']='DPSCF300028'
all_data_Z_FST_merge.loc[(all_data_Z_FST_merge['CHROM']=='DPSCF300028-2'), 'CHROM']='DPSCF300028'
all_data_Z_FST_merge.loc[(all_data_Z_FST_merge['CHROM']=='DPSCF300028-3'), 'CHROM']='DPSCF300028'
all_data_Z_FST_merge.loc[(all_data_Z_FST_merge['CHROM']=='DPSCF300044-1'), 'CHROM']='DPSCF300044'
all_data_Z_FST_merge.loc[(all_data_Z_FST_merge['CHROM']=='DPSCF300044-2'), 'CHROM']='DPSCF300044'


####################Genes in Z fst


all_genes_Z_FST=genes_in_outliers(all_data_Z_FST_merge)
#all_genes_Z_FST=pandas.merge(all_genes_Z_FST, Chrom, on=['CHROM'], how='inner')

genes_Z_FST = genes_in_outliers(data_Z_FST)
genes_Z_FST=pandas.merge(genes_Z_FST, Chrom, on=['CHROM'], how='inner')


Fst_Dxy_windows['outliers']=0
for window_index, window_1 in all_genes_Z_FST.iterrows():
	value=Fst_Dxy_windows.loc[(Fst_Dxy_windows['CHROM']==window_1.CHROM)&(Fst_Dxy_windows['BIN_START']==float(window_1.BIN_START))]['outliers']
	Fst_Dxy_windows.loc[(Fst_Dxy_windows['CHROM']==window_1.CHROM)&(Fst_Dxy_windows['BIN_START']==float(window_1.BIN_START)), 'outliers']=value+1


for i in [all_a_Z_FST, all_N_z_Z_FST, all_z_Z_FST]:
		mean=round(i[['Scaled_FST']].mean(), 4)
		std=round(i[['Scaled_FST']].std(), 4)
		print str(mean) + ' ± '+ str(std)

####################plots
head_stuff=['scaffold','source','feature','start','end','.','orentation','..','infor']
annotation=pandas.read_table("/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_genes.gff",skiprows=0, names=head_stuff)
annotation=annotation[['scaffold','source','feature','start','end','infor','orentation']]
annotation_gene=annotation[annotation.feature=='gene']
annotation_gene['ID']=annotation_gene.infor.str.split(';').str.get(0).str.split('=').str.get(1)

annotation_gene=annotation_gene[['scaffold', 'start', 'end', 'orentation','ID']]
annotation_gene=annotation_gene.rename(columns={'scaffold': 'CHROM', 'start': 'BIN_START', 'end': 'BIN_END'})
annotation_gene['Sites']=0
annotation_gene['Sites']=annotation_gene['BIN_END']-annotation_gene['BIN_START']
annotation_gene=annotation_gene[['CHROM','BIN_START','BIN_END','Sites','ID']]
	
def callcualte_gene_FST_from_site_FST(population_pi, dataframe, sites_list, output):
	win_all_sites=dataframe[['CHROM','BIN_START','BIN_END',sites_list]]
	win_sites=population_pi
	df_merge=pandas.DataFrame(columns=['CHROM','POS','WEIR_AND_COCKERHAM_FST','BIN_START','BIN_END','Sites', 'ID'])
	l=list(set(win_all_sites.CHROM))
	l.sort()
	for i in l:
		chrom_set=win_all_sites[win_all_sites['CHROM']==i]
		chrom_sites=win_sites[win_sites['CHROM']==i]
		new = pandas.merge(chrom_set, chrom_sites, on=['CHROM'],how='inner')
		new = new.query('BIN_START <= POS and BIN_END >= POS')
		df_merge=pandas.concat([new, df_merge], ignore_index=True)
	out=pandas.DataFrame({output :df_merge.groupby(['CHROM','BIN_START','BIN_END',sites_list], as_index=True)['WEIR_AND_COCKERHAM_FST'].sum()}).reset_index()
	out[output]=out[output]/out[sites_list]
	return out


def callcualte_windowed_pi_from_site_pi(population_pi, dataframe, sites_list, output):
	win_all_sites=dataframe[['CHROM','BIN_START','BIN_END',sites_list,'ID']]
	win_sites=population_pi
	df_merge=pandas.DataFrame(columns=['CHROM','POS','PI','BIN_START','BIN_END','Sites', 'ID'])
	l=list(set(win_all_sites.CHROM))
	l.sort()
	for i in l:
		chrom_set=win_all_sites[win_all_sites['CHROM']==i]
		chrom_sites=win_sites[win_sites['CHROM']==i]
		new = pandas.merge(chrom_set, chrom_sites, on=['CHROM'],how='inner')
		new = new.query('BIN_START <= POS and BIN_END >= POS')
		df_merge=pandas.concat([new, df_merge], ignore_index=True)
	out=pandas.DataFrame({output+'_gene_pi' :df_merge.groupby(['CHROM','BIN_START','BIN_END',sites_list], as_index=True)['PI'].sum()}).reset_index()
	out[output+'_gene_pi']=out[output+'_gene_pi']/out[sites_list]
	return out



#######################pop gen stats for genes
#Fst_sites=pandas.read_table('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/Fst_window_site.weir.fst')
#Fst_sites['WEIR_AND_COCKERHAM_FST']=Fst_sites['WEIR_AND_COCKERHAM_FST'].clip(lower=0)
#Fst_sites[Fst_sites['WEIR_AND_COCKERHAM_FST']  > Fst_sites.WEIR_AND_COCKERHAM_FST.dropna().quantile(0.99)]
#annotation_gene_Fst_sites=callcualte_gene_FST_from_site_FST(Fst_sites, annotation_gene, 'Sites', 'gene_MEAN_FST')
#annotation_gene_Fst_sites=pandas.merge(annotation_gene_Fst_sites, annotation_gene, on=['CHROM', 'BIN_START', 'BIN_END', 'Sites'], how='inner')
#annotation_gene_Fst_sites=pandas.merge(annotation_gene_Fst_sites, annotation_gene, on=['CHROM', 'BIN_START', 'BIN_END', 'Sites'], how='inner')

#East ="/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/East.PI.sites.sites.pi"
#West ="/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/West.PI.sites.sites.pi"
#East_site_PI_table =pandas.read_table(East)
#West_site_PI_table =pandas.read_table(West)


##calculate gene wise pi for east and west
#annotation_East_site_PI=callcualte_windowed_pi_from_site_pi(East_site_PI_table, annotation_gene, 'Sites', 'East')
#annotation_West_site_PI=callcualte_windowed_pi_from_site_pi(West_site_PI_table, annotation_gene, 'Sites', 'West')


##merge gene stats
annotation_gene_Fst_sites=pandas.read_csv('Fst_genes_DB.csv')
annotation_gene_Fst_sites=pandas.merge(annotation_gene_Fst_sites[["CHROM","BIN_START","BIN_END","Sites","gene_MEAN_FST"]], annotation_gene, on=['CHROM','BIN_START','BIN_END','Sites'], how='outer')
#annotation_gene_Fst_sites=pandas.merge(annotation_gene_Fst_sites, annotation_West_site_PI, on=['CHROM','BIN_START','BIN_END','Sites'], how='outer')
#annotation_gene_Fst_sites=pandas.merge(annotation_gene_Fst_sites, annotation_East_site_PI, on=['CHROM','BIN_START','BIN_END','Sites'], how='outer')
#


genes_from_paper=[ 'DPOGS200868', 'DPOGS201379', 'DPOGS211203', 'DPOGS202675', 'DPOGS215054', 'DPOGS211604']

annotation_gene_Fst_sites[annotation_gene_Fst_sites.ID.isin(genes_from_paper)]

annotation_gene_Fst_sites.loc[(annotation_gene_Fst_sites['CHROM']=='DPSCF300001')&(annotation_gene_Fst_sites['BIN_START']<5820000), 'CHROM']='DPSCF300001-1'
annotation_gene_Fst_sites.loc[(annotation_gene_Fst_sites['CHROM']=='DPSCF300001')&(annotation_gene_Fst_sites['BIN_START']>5820000), 'CHROM']='DPSCF300001-2'
annotation_gene_Fst_sites.loc[(annotation_gene_Fst_sites['CHROM']=='DPSCF300028')&(annotation_gene_Fst_sites['BIN_START']<760000), 'CHROM']='DPSCF300028-1'
annotation_gene_Fst_sites.loc[(annotation_gene_Fst_sites['CHROM']=='DPSCF300028')&(annotation_gene_Fst_sites['BIN_START']>1805000), 'CHROM']='DPSCF300028-3'
annotation_gene_Fst_sites.loc[(annotation_gene_Fst_sites['CHROM']=='DPSCF300028')&(annotation_gene_Fst_sites['BIN_START']>760000)&(annotation_gene_Fst_sites['BIN_START']<1805000), 'CHROM']='DPSCF300028-2'
annotation_gene_Fst_sites.loc[(annotation_gene_Fst_sites['CHROM']=='DPSCF300044')&(annotation_gene_Fst_sites['BIN_START']<290000), 'CHROM']='DPSCF300044-1'
annotation_gene_Fst_sites.loc[(annotation_gene_Fst_sites['CHROM']=='DPSCF300044')&(annotation_gene_Fst_sites['BIN_START']>290000), 'CHROM']='DPSCF300044-2'

annotation_gene_Fst_sites.loc[(annotation_gene_Fst_sites['CHROM']=='DPSCF300001-1')&(annotation_gene_Fst_sites['BIN_START']>3827349)&(annotation_gene_Fst_sites['BIN_END']<3882816), 'CHROM']='DPSCF300001-Z_neo'


Gene_fst_windows=Fst_Dxy_windows[['CHROM', 'Auto_Allo', 'Chr_number']]
annotation_gene_Fst_sites_2=pandas.merge(annotation_gene_Fst_sites[["CHROM","BIN_START","BIN_END","Sites","gene_MEAN_FST", 'ID']], Gene_fst_windows, on=['CHROM'], how='outer')
annotation_gene_Fst_sites_2=annotation_gene_Fst_sites_2.drop_duplicates(keep='first')
############outliers


z=annotation_gene_Fst_sites_2[(annotation_gene_Fst_sites_2['Auto_Allo'] == 'Z_anc')]
N_z=annotation_gene_Fst_sites_2[(annotation_gene_Fst_sites_2['Auto_Allo'] == 'Z_neo')]
a=annotation_gene_Fst_sites_2[(annotation_gene_Fst_sites_2['Auto_Allo'] == 'A')]
gene_a_gene_MEAN_FST =a[(a['gene_MEAN_FST']  > a.gene_MEAN_FST.dropna().quantile(0.99))]
gene_N_z_gene_MEAN_FST =N_z[(N_z['gene_MEAN_FST']  > N_z.gene_MEAN_FST.dropna().quantile(0.99))]
gene_z_gene_MEAN_FST =z[(z['gene_MEAN_FST']  > z.gene_MEAN_FST.dropna().quantile(0.99))]

##############Tajmas_D_outliers
z=Fst_Dxy_windows[(Fst_Dxy_windows['Auto_Allo'] == 'Z_anc')]
N_z=Fst_Dxy_windows[(Fst_Dxy_windows['Auto_Allo'] == 'Z_neo')]
a=Fst_Dxy_windows[(Fst_Dxy_windows['Auto_Allo'] == 'A')]


z['East_TajimaD_Z']=(z['East_TajimaD']-np.nanmean(z['East_TajimaD']))/np.nanstd(z['East_TajimaD'])
N_z['East_TajimaD_Z']=(N_z['East_TajimaD']-np.nanmean(N_z['East_TajimaD']))/np.nanstd(N_z['East_TajimaD'])
a['East_TajimaD_Z']=(a['East_TajimaD']-np.nanmean(a['East_TajimaD']))/np.nanstd(a['East_TajimaD'])


z['West_TajimaD_Z']=(z['West_TajimaD']-np.nanmean(z['West_TajimaD']))/np.nanstd(z['West_TajimaD'])
N_z['West_TajimaD_Z']=(N_z['West_TajimaD']-np.nanmean(N_z['West_TajimaD']))/np.nanstd(N_z['West_TajimaD'])
a['West_TajimaD_Z']=(a['West_TajimaD']-np.nanmean(a['West_TajimaD']))/np.nanstd(a['West_TajimaD'])


East_a_Td =a[(a['East_TajimaD_Z']  < a.East_TajimaD_Z.dropna().quantile(0.05))]
East_N_z_Td =N_z[(N_z['East_TajimaD_Z']  < N_z.East_TajimaD_Z.dropna().quantile(0.05))]
East_z_Td =z[(z['East_TajimaD_Z']  < z.East_TajimaD.dropna().quantile(0.05))]

West_a_Td =a[(a['West_TajimaD_Z']  < a.West_TajimaD_Z.dropna().quantile(0.05))]
West_N_z_Td =N_z[(N_z['West_TajimaD_Z']  < N_z.West_TajimaD_Z.dropna().quantile(0.05))]
West_z_Td =z[(z['West_TajimaD_Z']  < z.West_TajimaD_Z.dropna().quantile(0.05))]

East_Td =pandas.concat([East_a_Td ,East_N_z_Td, East_z_Td], ignore_index=True)
West_Td =pandas.concat([West_a_Td ,West_N_z_Td, West_z_Td], ignore_index=True)


East_Td.loc[(East_Td['CHROM']=='DPSCF300044-1'), 'CHROM']='DPSCF300044'
East_Td.loc[(East_Td['CHROM']=='DPSCF300001-1'), 'CHROM']='DPSCF300001'
East_Td.loc[(East_Td['CHROM']=='DPSCF300028-1'), 'CHROM']='DPSCF300028'
East_Td.loc[(East_Td['CHROM']=='DPSCF300001-2'), 'CHROM']='DPSCF300001'


West_Td.loc[(West_Td['CHROM']=='DPSCF300044-1'), 'CHROM']='DPSCF300044'
West_Td.loc[(West_Td['CHROM']=='DPSCF300001-1'), 'CHROM']='DPSCF300001'
West_Td.loc[(West_Td['CHROM']=='DPSCF300028-1'), 'CHROM']='DPSCF300028'
West_Td.loc[(West_Td['CHROM']=='DPSCF300001-2'), 'CHROM']='DPSCF300001'

##############outlier genes
Fst_genes_outliers=pandas.read_csv('gene_Fst_outliers.csv', header=None,names=['gene'])
Td_East_genes=genes_in_outliers(East_Td)
Td_West_genes=genes_in_outliers(West_Td)


common_Fst_East=pandas.merge(Fst_genes_outliers, Td_East_genes, on=['gene'], how='inner')
common_Fst_West=pandas.merge(Fst_genes_outliers, Td_West_genes, on=['gene'], how='inner')

#######################

sns.palplot(sns.color_palette("colorblind", 12))


manhattan_plot_fst_dxy_pi(Fst_Dxy_windows, ['WEIR_AND_COCKERHAM_MEAN_FST'], ['East_West_dxy'], ['East_final_pi','West_final_pi'], ['East_TajimaD','West_TajimaD'],['outliers'],CHROM_dict, '#050505')

#manhattan_plot_fst_dxy_pi(sites_filter_0, fst, dxy,pi,taj,fixed,chromosomes_dict,colour)
manhattan_plot_fst_dxy_pi(Fst_Dxy_windows, ['WEIR_AND_COCKERHAM_MEAN_FST'], ['East_West_dxy'], ['East_final_pi','West_final_pi'], ['East_TajimaD','West_TajimaD'],['East_West_fixed','East_West_private_a','East_West_private_b','East_West_shared'],CHROM_dict, '#050505')


plt.savefig('Fst_Dxy_windows.png')
sns.boxplot(data=FST_windows[['East_TajimaD','West_TajimaD']],showfliers=False, palette=sns.color_palette(['#d90d0d','#606060'])).set(ylim=(-3.0, 3.0))

sns.boxplot(data=Fst_Dxy_windows[['East_final_pi','West_final_pi']],showfliers=False, palette=sns.color_palette(['#d90d0d','#606060']))

sns.jointplot('West_final_pi','East_final_pi', data=Fst_Dxy_windows, kind="reg")

sns.jointplot('WEIR_AND_COCKERHAM_MEAN_FST','East_West_dxy', data=Fst_Dxy_windows, kind="reg")

sns.distplot(FST_windows[['East_TajimaD','West_TajimaD']], hist=False, kde_kws={"shade": True})

#test for signisficance mannwhitneyu test
sns.kdeplot(Fst_Dxy_windows.East_TajimaD, shade=True) 
sns.kdeplot(Fst_Dxy_windows.West_TajimaD, shade=True) 

sns.jointplot('East_TajimaD', 'West_TajimaD', data=Fst_Dxy_windows, kind="kde", space=0)
sns.jointplot('East_TajimaD', 'West_TajimaD', data=Fst_Dxy_windows, kind="reg", space=0)

sns.jointplot('East_West_dxy', 'Avg_pi', data=Fst_Dxy_windows, kind="reg", space=0)


scipy.stats.mannwhitneyu(list(FST_windows.East_TajimaD), list(FST_windows.West_TajimaD))
#MannwhitneyuResult(statistic=106843526.0, pvalue=0.0)

sns.kdeplot(z['Z_FST'], shade=True) 

corr_plot(Fst_Dxy_windows[['WEIR_AND_COCKERHAM_MEAN_FST','Carpinteria_VS_Oceano_MEAN_FST','Big_Sur_VS_Oceano_MEAN_FST','Big_Sur_VS_Carpinteria_MEAN_FST']])
corr_pvalue(Fst_Dxy_windows[['WEIR_AND_COCKERHAM_MEAN_FST','Carpinteria_VS_Oceano_MEAN_FST','Big_Sur_VS_Oceano_MEAN_FST','Big_Sur_VS_Carpinteria_MEAN_FST']])

##################


MK_test=pandas.read_csv("/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/FST_outliers_GO.csv")


x = list(MK_test['classicFisher'])
y = list(MK_test['Term'])
z = list(MK_test['Significant']*4)

c = list(MK_test['Significant']/MK_test['Annotated'])

sc=plt.scatter(x, y, s=z, c=c, cmap="viridis", edgecolors="grey", linewidth=0.01)
plt.colorbar(sc)
plt.show()


0.0079 ± 0.0027
