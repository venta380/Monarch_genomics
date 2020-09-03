import sys
import string
import numpy as np
import fileinput
import pandas
import personal_popgen
import argparse
import os
import functools
import subprocess


pwd='/scratch/vt20265/bam_files/'
os.chdir(pwd)

files=subprocess.check_output("ls *_N_sites", shell=True, universal_newlines=True).split()


population_names=set([i.split('_')[0] for i in files])
type_names=set(['all' if 'intergenic' not in i.split('_') for i in files ])




for i in population_names:
	files_pop=set([j  for j in files if i in j.split('_')])
	files_pop_df=pandas.DataFrame()
	for n in files_pop: 
		size=n.split('_')[1]
		new=pandas.read_csv(n)
		files_pop_df=files_pop_df.append(new, ignore_index=True)
	files_pop_df.to_csv(i+'_'+size+'_N_sites', index=False, sep='\t', header=False, )

