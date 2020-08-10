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
import random 
from itertools import izip



fasta_DPlex=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_scaffolds.fa')

fasta_stm163_1=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/stm163_1.fa')
fasta_stm146_1=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/stm146_1.fa')
fasta_T9_1=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/T9_1.fa')
fasta_T14_1=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/T14_1.fa')
fasta_NJ203_1=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/NJ203_1.fa')
fasta_NJ116_1=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/NJ116_1.fa')
fasta_NJ1_1=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/NJ1_1.fa')
fasta_HI023_1=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/HI023_1.fa')
fasta_HI033_1=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/HI033_1.fa')
fasta_mex986_1=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/mex986_1.fa')
fasta_mex919_1=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/mex919_1.fa')
fasta_mex915_1=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/mex915_1.fa')
fasta_mex536_1=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/mex536_1.fa')
fasta_mex1527_1=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/mex1527_1.fa')

fasta_stm163_2=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/stm163_2.fa')
fasta_stm146_2=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/stm146_2.fa')
fasta_T9_2=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/T9_2.fa')
fasta_T14_2=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/T14_2.fa')
fasta_NJ203_2=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/NJ203_2.fa')
fasta_NJ116_2=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/NJ116_2.fa')
fasta_NJ1_2=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/NJ1_2.fa')
fasta_HI023_2=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/HI023_2.fa')
fasta_HI033_2=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/HI033_2.fa')
fasta_mex986_2=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/mex986_2.fa')
fasta_mex919_2=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/mex919_2.fa')
fasta_mex915_2=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/mex915_2.fa')
fasta_mex536_2=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/mex536_2.fa')
fasta_mex1527_2=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/mex1527_2.fa')



for i in fasta_DPlex.keys():
	start_list=range(0,len(fasta_DPlex[i]), 50000)
	for start in start_list:
		end=start+50000
		if len(str(fasta_DPlex[i][start:end])) > 40000:
			ns=fasta_DPlex[i][start:end].count('N')
			gs=fasta_DPlex[i][start:end].count('G')
			cs=fasta_DPlex[i][start:end].count('C')
			length=len(fasta_DPlex[i][start:end])
			GC_cont=(float(gs)+float(cs))/(length-ns)
			f = open('/crex/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/windows/'+str(i)+"_"+str(start)+"_"+str(GC_cont)+'.fa', "w")
			f.write(">stm163_1"+"\n")
			f.write(str(fasta_stm163_1[i][start:end])+"\n")
			f.write(">stm146_1"+"\n")
			f.write(str(fasta_stm146_1[i][start:end])+"\n")
			f.write(">T9_1"+"\n")
			f.write(str(fasta_T9_1[i][start:end])+"\n")
			f.write(">T14_1"+"\n")
			f.write(str(fasta_T14_1[i][start:end])+"\n")
			f.write(">NJ203_1"+"\n")
			f.write(str(fasta_NJ203_1[i][start:end])+"\n")
			f.write(">NJ116_1"+"\n")
			f.write(str(fasta_NJ116_1[i][start:end])+"\n")
			f.write(">NJ1_1"+"\n")
			f.write(str(fasta_NJ1_1[i][start:end])+"\n")
			f.write(">HI023_1"+"\n")
			f.write(str(fasta_HI023_1[i][start:end])+"\n")
			f.write(">HI033_1"+"\n")
			f.write(str(fasta_HI033_1[i][start:end])+"\n")
			f.write(">mex986_1"+"\n")
			f.write(str(fasta_mex986_1[i][start:end])+"\n")
			f.write(">mex919_1"+"\n")
			f.write(str(fasta_mex919_1[i][start:end])+"\n")
			f.write(">mex915_1"+"\n")
			f.write(str(fasta_mex915_1[i][start:end])+"\n")
			f.write(">mex536_1"+"\n")
			f.write(str(fasta_mex536_1[i][start:end])+"\n")
			f.write(">mex1527_1"+"\n")
			f.write(str(fasta_mex1527_1[i][start:end])+"\n")
			f.write(">stm163_2"+"\n")
			f.write(str(fasta_stm163_2[i][start:end])+"\n")
			f.write(">stm146_2"+"\n")
			f.write(str(fasta_stm146_2[i][start:end])+"\n")
			f.write(">T9_2"+"\n")
			f.write(str(fasta_T9_2[i][start:end])+"\n")
			f.write(">T14_2"+"\n")
			f.write(str(fasta_T14_2[i][start:end])+"\n")
			f.write(">NJ203_2"+"\n")
			f.write(str(fasta_NJ203_2[i][start:end])+"\n")
			f.write(">NJ116_2"+"\n")
			f.write(str(fasta_NJ116_2[i][start:end])+"\n")
			f.write(">NJ1_2"+"\n")
			f.write(str(fasta_NJ1_2[i][start:end])+"\n")
			f.write(">HI023_2"+"\n")
			f.write(str(fasta_HI023_2[i][start:end])+"\n")
			f.write(">HI033_2"+"\n")
			f.write(str(fasta_HI033_2[i][start:end])+"\n")
			f.write(">mex986_2"+"\n")
			f.write(str(fasta_mex986_2[i][start:end])+"\n")
			f.write(">mex919_2"+"\n")
			f.write(str(fasta_mex919_2[i][start:end])+"\n")
			f.write(">mex915_2"+"\n")
			f.write(str(fasta_mex915_2[i][start:end])+"\n")
			f.write(">mex536_2"+"\n")
			f.write(str(fasta_mex536_2[i][start:end])+"\n")
			f.write(">mex1527_2"+"\n")
			f.write(str(fasta_mex1527_2[i][start:end])+"\n")
			f.close()

