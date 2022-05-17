import sys
import os
import numpy
import dadi
import pylab
from datetime import datetime
import Optimize_Functions
import Models_2D
import os
import numpy
import dadi
import Plotting_Functions
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum
import pickle


def bottlegrowth(params, ns, pts):
    """
    params = (nuB,nuF,T)
    ns = (n1,n2)
    Instantanous size change followed by exponential growth with no population
    split.
    nuB: Ratio of population size after instantanous change to ancient
         population size
    nuF: Ratio of contempoary to ancient population size
    T: Time in the past at which instantaneous change happened and growth began
       (in units of 2*Na generations) 
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nuB,nuF,T = params
    return bottlegrowth_split_mig((nuB,nuF,0,T,0), ns, pts)

def bottlegrowth_split(params, ns, pts):
    """
    params = (nuB,nuF,T,Ts)
    ns = (n1,n2)
    Instantanous size change followed by exponential growth then split.
    nuB: Ratio of population size after instantanous change to ancient
         population size
    nuF: Ratio of contempoary to ancient population size
    T: Time in the past at which instantaneous change happened and growth began
       (in units of 2*Na generations) 
    Ts: Time in the past at which the two populations split.
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nuB,nuF,T,Ts = params
    return bottlegrowth_split_mig((nuB,nuF,0,T,Ts), ns, pts)

def bottlegrowth_split_mig(params, ns, pts):
    """
    params = (nuB,nuF,m,T,Ts)
    ns = (n1,n2)
    Instantanous size change followed by exponential growth then split with
    migration.
    nuB: Ratio of population size after instantanous change to ancient
         population size
    nuF: Ratio of contempoary to ancient population size
    m: Migration rate between the two populations (2*Na*m).
    T: Time in the past at which instantaneous change happened and growth began
       (in units of 2*Na generations) 
    Ts: Time in the past at which the two populations split.
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nuB,nuF,m,T,Ts = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)

    nu_func = lambda t: nuB*numpy.exp(numpy.log(nuF/nuB) * t/T)
    phi = Integration.one_pop(phi, xx, T-Ts, nu_func)

    phi = PhiManip.phi_1D_to_2D(xx, phi)
    nu0 = nu_func(T-Ts)
    nu_func = lambda t: nu0*numpy.exp(numpy.log(nuF/nu0) * t/Ts)
    phi = Integration.two_pops(phi, xx, Ts, nu_func, nu_func, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def split_mig(params, ns, pts):
    """
    params = (nu1,nu2,T,m)
    ns = (n1,n2)
    Split into two populations of specifed size, with migration.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    m: Migration rate between populations (2*Na*m)
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1,nu2,T,m = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


#
#dd = dadi.Misc.make_data_dict_vcf("/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/Autosome.recode.vcf","/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/pop_All.txt")
#fs = dadi.Spectrum.from_data_dict( dd , [ "ALL"], projections =[43] ,polarized = False )
#file = open('Single_pop.txt', 'wb')
#pickle.dump(fs, file)


sfs='/home/venkat/bin/dadi_pipeline/Two_Population_Pipeline_Auto/Single_pop.txt'
file_sfs= open(sfs, 'r')
fs=pickle.load(file_sfs)

proj = [43]


pts = [50,60,70]

rounds = 4

reps = [10,20,30,40]
maxiters = [3,5,10,15]
folds = [3,2,2,1]

fs_folded = True


import Demographics1D


prefix = "ALL"

#snm
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "snm", Demographics1D.snm, rounds, 1, fs_folded=fs_folded,
                                        optimizer="log",reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1")

#two_epoch
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "two_epoch", Demographics1D.two_epoch, rounds, 2, fs_folded=fs_folded,
                                        optimizer="log",reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu,T")

#growth
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "growth", Demographics1D.growth, rounds, 2, fs_folded=fs_folded,
                                        optimizer="log",reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu,T")

#bottlegrowth
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "bottlegrowth", Demographics1D.bottlegrowth, rounds, 3, fs_folded=fs_folded,
                                        optimizer="log",reps=reps, maxiters=maxiters, folds=folds, param_labels = "nuB,nuF,T")

##bottlegrowth
#Optimize_Functions.Optimize_Routine(fs, pts, prefix, "bottlegrowth", Demographics1D.bottlegrowth, rounds, 3, fs_folded=fs_folded,
#                                        optimizer="log",reps=reps, maxiters=maxiters, folds=folds, param_labels = "nuB,nuF,T")


#three_epoch
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "three_epoch", Demographics1D.three_epoch, rounds, 4, fs_folded=fs_folded,
                                        optimizer="log",reps=reps, maxiters=maxiters, folds=folds, param_labels = "nuB,nuF,TB,TF")



import Demographics2D




sfs='/home/venkat/bin/dadi_pipeline/Two_Population_Pipeline_Auto/sfs_Auto.txt'
file_sfs= open(sfs, 'r')
fs=pickle.load(file_sfs)
proj = [29,14]
pts = [50,60,70]
rounds = 4
reps = [10,20,30,40]
maxiters = [3,5,10,15]
folds = [3,2,2,1]
fs_folded = True


prefix="Two_pop_"

#bottlegrowth
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "bottlegrowth", Demographics2D.bottlegrowth, rounds, 3, fs_folded=fs_folded,
                                        optimizer="log", reps=reps, maxiters=maxiters, folds=folds, param_labels = "nuB,nuF,T")

#bottlegrowth_split
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "bottlegrowth_split", Demographics2D.bottlegrowth_split, rounds, 4, fs_folded=fs_folded,
                                        optimizer="log", reps=reps, maxiters=maxiters, folds=folds, param_labels = "nuB,nuF,T,Ts")


#bottlegrowth_split_mig
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "bottlegrowth_split_mig", Demographics2D.bottlegrowth_split_mig, rounds, 5, fs_folded=fs_folded,
                                        optimizer="log", reps=reps, maxiters=maxiters, folds=folds, param_labels = "nuB,nuF,m,T,Ts")


#split_mig
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "split_mig", Demographics2D.split_mig, rounds, 4, fs_folded=fs_folded,
                                        optimizer="log", reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1,nu2,T,m")


#IM
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "IM", Demographics2D.IM, rounds, 6, fs_folded=fs_folded,
                                        optimizer="log", reps=reps, maxiters=maxiters, folds=folds, param_labels = "s,nu1,nu2,T,m12,m21")


#IM_pre
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "IM_pre", Demographics2D.IM_pre, rounds, 8, fs_folded=fs_folded,
                                        optimizer="log", reps=reps, maxiters=maxiters, folds=folds, param_labels = "nuPre,TPre,s,nu1,nu2,T,m12,m21")



emp_params = [1.7745,6.3606,10.9812,0.9189,0.2505]

model_fit = Plotting_Functions.Fit_Empirical(fs, pts, prefix, "bottlegrowth_split_mig", bottlegrowth_split_mig, emp_params, fs_folded=True)

Plotting_Functions.Plot_2D(fs, model_fit, prefix, "bottlegrowth_split_mig")



bottlegrowth_split_mig  Round_4_Replicate_26    -3306.97        6623.94 4637.08 980321.63       1.7745,6.3606,10.9812,0.9189,0.2505


L=11288846.39611162
theta=980321.63
u=2.9*10e-9
Ne=theta/(4*u*L)


(nuB,nuF,m,T,Ts)

nuB=1.7745*Ne
nuF=6.3606*Ne
m=10.9812/(2*Ne)
T=0.9189*(2*Ne)*0.3
Ts=0.2505*(2*Ne)*0.3




