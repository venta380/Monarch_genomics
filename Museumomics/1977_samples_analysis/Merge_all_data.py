import sys
import os
import string
import pandas
import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns
import itertools
import math
import time
import sys
import personal_popgen
import re
import random 
from Bio import SeqIO
from Bio.Seq import Seq
import argparse

 

West_all_df=pandas.DataFrame()
West_intergenic_df=pandas.DataFrame()
West_codon1_df=pandas.DataFrame()
West_codon2_df=pandas.DataFrame()
West_codon3_df=pandas.DataFrame()
West_4fold_df=pandas.DataFrame()
East_all_df=pandas.DataFrame()
East_intergenic_df=pandas.DataFrame()
East_codon1_df=pandas.DataFrame()
East_codon2_df=pandas.DataFrame()
East_codon3_df=pandas.DataFrame()
East_4fold_df=pandas.DataFrame()
_1977_all_df=pandas.DataFrame()
_1977_intergenic_df=pandas.DataFrame()
_1977_codon1_df=pandas.DataFrame()
_1977_codon2_df=pandas.DataFrame()
_1977_codon3_df=pandas.DataFrame()
_1977_4fold_df=pandas.DataFrame()

for i in range(0,20):
	West_all_df_name="West_all_"+str(i)+"_stats.csv"
	West_intergenic_df_name="West_intergenic_"+str(i)+"_stats.csv"
	West_codon1_df_name="West_codon1_"+str(i)+"_stats.csv"
	West_codon2_df_name="West_codon2_"+str(i)+"_stats.csv"
	West_codon3_df_name="West_codon3_"+str(i)+"_stats.csv"
	West_4fold_df_name="West_4fold_"+str(i)+"_stats.csv"
	East_all_df_name="East_all_"+str(i)+"_stats.csv"
	East_intergenic_df_name="East_intergenic_"+str(i)+"_stats.csv"
	East_codon1_df_name="East_codon1_"+str(i)+"_stats.csv"
	East_codon2_df_name="East_codon2_"+str(i)+"_stats.csv"
	East_codon3_df_name="East_codon3_"+str(i)+"_stats.csv"
	East_4fold_df_name="East_4fold_"+str(i)+"_stats.csv"
	_1977_all_df_name="1977_all_"+str(i)+"_stats.csv"
	_1977_intergenic_df_name="1977_intergenic_"+str(i)+"_stats.csv"
	_1977_codon1_df_name="1977_codon1_"+str(i)+"_stats.csv"
	_1977_codon2_df_name="1977_codon2_"+str(i)+"_stats.csv"
	_1977_codon3_df_name="1977_codon3_"+str(i)+"_stats.csv"
	_1977_4fold_df_name="1977_4fold_"+str(i)+"_stats.csv"
	West_all_df=pandas.concat([West_all_df,pandas.read_csv(West_all_df_name, skiprows=1,names=['CHROM','BIN_START', 'West_all_Pi','West_all_Tw','West_all_Td', 'West_all_Nsites'],header=None)])
	West_intergenic_df=pandas.concat([West_intergenic_df,pandas.read_csv(West_intergenic_df_name, skiprows=1,names=['CHROM','BIN_START', 'West_intergenic_Pi','West_intergenic_Tw','West_intergenic_Td', 'West_intergenic_Nsites'],header=None)])
	West_codon1_df=pandas.concat([West_codon1_df,pandas.read_csv(West_codon1_df_name, skiprows=1,names=['CHROM','BIN_START', 'West_codon1_Pi','West_codon1_Tw','West_codon1_Td', 'West_codon1_Nsites'],header=None)])
	West_codon2_df=pandas.concat([West_codon2_df,pandas.read_csv(West_codon2_df_name, skiprows=1,names=['CHROM','BIN_START', 'West_codon2_Pi','West_codon2_Tw','West_codon2_Td', 'West_codon2_Nsites'],header=None)])
	West_codon3_df=pandas.concat([West_codon3_df,pandas.read_csv(West_codon3_df_name, skiprows=1,names=['CHROM','BIN_START', 'West_codon3_Pi','West_codon3_Tw','West_codon3_Td', 'West_codon3_Nsites'],header=None)])
	West_4fold_df=pandas.concat([West_4fold_df,pandas.read_csv(West_4fold_df_name, skiprows=1,names=['CHROM','BIN_START', 'West_4fold_Pi','West_4fold_Tw','West_4fold_Td', 'West_4fold_Nsites'],header=None)])
	East_all_df=pandas.concat([East_all_df,pandas.read_csv(East_all_df_name, skiprows=1,names=['CHROM','BIN_START', 'East_all_Pi','East_all_Tw','East_all_Td', 'East_all_Nsites'],header=None)])
	East_intergenic_df=pandas.concat([East_intergenic_df,pandas.read_csv(East_intergenic_df_name, skiprows=1,names=['CHROM','BIN_START', 'East_intergenic_Pi','East_intergenic_Tw','East_intergenic_Td', 'East_intergenic_Nsites'],header=None)])
	East_codon1_df=pandas.concat([East_codon1_df,pandas.read_csv(East_codon1_df_name, skiprows=1,names=['CHROM','BIN_START', 'East_codon1_Pi','East_codon1_Tw','East_codon1_Td', 'East_codon1_Nsites'],header=None)])
	East_codon2_df=pandas.concat([East_codon2_df,pandas.read_csv(East_codon2_df_name, skiprows=1,names=['CHROM','BIN_START', 'East_codon2_Pi','East_codon2_Tw','East_codon2_Td', 'East_codon2_Nsites'],header=None)])
	East_codon3_df=pandas.concat([East_codon3_df,pandas.read_csv(East_codon3_df_name, skiprows=1,names=['CHROM','BIN_START', 'East_codon3_Pi','East_codon3_Tw','East_codon3_Td', 'East_codon3_Nsites'],header=None)])
	East_4fold_df=pandas.concat([East_4fold_df,pandas.read_csv(East_4fold_df_name, skiprows=1,names=['CHROM','BIN_START', 'East_4fold_Pi','East_4fold_Tw','East_4fold_Td', 'East_4fold_Nsites'],header=None)])
	_1977_all_df=pandas.concat([_1977_all_df,pandas.read_csv(_1977_all_df_name, skiprows=1,names=['CHROM','BIN_START', '1977_all_Pi','1977_all_Tw','1977_all_Td', '1977_all_Nsites'],header=None)])
	_1977_intergenic_df=pandas.concat([_1977_intergenic_df,pandas.read_csv(_1977_intergenic_df_name, skiprows=1,names=['CHROM','BIN_START', '1977_intergenic_Pi','1977_intergenic_Tw','1977_intergenic_Td', '1977_intergenic_Nsites'],header=None)])
	_1977_codon1_df=pandas.concat([_1977_codon1_df,pandas.read_csv(_1977_codon1_df_name, skiprows=1,names=['CHROM','BIN_START', '1977_codon1_Pi','1977_codon1_Tw','1977_codon1_Td', '1977_codon1_Nsites'],header=None)])
	_1977_codon2_df=pandas.concat([_1977_codon2_df,pandas.read_csv(_1977_codon2_df_name, skiprows=1,names=['CHROM','BIN_START', '1977_codon2_Pi','1977_codon2_Tw','1977_codon2_Td', '1977_codon2_Nsites'],header=None)])
	_1977_codon3_df=pandas.concat([_1977_codon3_df,pandas.read_csv(_1977_codon3_df_name, skiprows=1,names=['CHROM','BIN_START', '1977_codon3_Pi','1977_codon3_Tw','1977_codon3_Td', '1977_codon3_Nsites'],header=None)])
	_1977_4fold_df=pandas.concat([_1977_4fold_df,pandas.read_csv(_1977_4fold_df_name, skiprows=1,names=['CHROM','BIN_START', '1977_4fold_Pi','1977_4fold_Tw','1977_4fold_Td', '1977_4fold_Nsites'],header=None)])




