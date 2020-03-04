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


def snm(notused, ns, pts):
    """
    ns = (n1,n2)
    Standard neutral model, populations never diverge.
    """
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def no_divergence(notused, ns, pts):
    """
    Standard neutral model, populations never diverge.
    """

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def no_mig(params, ns, pts):
    """
    Split into two populations, no migration.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations)
    """
    nu1, nu2, T = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def sym_mig(params, ns, pts):
    """
    Split into two populations, with symmetric migration.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations)
    m: Migration rate between populations (2*Na*m)
    """
    nu1, nu2, m, T = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def asym_mig(params, ns, pts):
    """
    Split into two populations, with different migration rates.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations)
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
        """
    nu1, nu2, m12, m21, T = params
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m12, m21=m21)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))

    return fs



#dd = dadi.Misc.make_data_dict_vcf("/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/good_pos.vcf.gz","/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/pop.txt")
#fs = dadi.Spectrum.from_data_dict( dd , [ "West"  , "East" ] ,projections =[29,14] ,polarized = False )
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


#no dovergence
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "no_divergence", Models_2D.no_divergence, rounds, 1, fs_folded=fs_folded,
                                        optimizer="log", maxiters=maxiters, folds=folds, param_labels = "nu1")


#Split into two populations, no migration.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "no_mig", Models_2D.no_mig, rounds, 3, fs_folded=fs_folded,
                                        optimizer="log", maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, T")


#Split into two populations, with continuous asymmetric migration.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "asym_mig", Models_2D.asym_mig, rounds, 5, fs_folded=fs_folded,
                                        optimizer="log", maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, T")

# Split into two populations, with continuous symmetric migration.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sym_mig", Models_2D.sym_mig, rounds, 4, fs_folded=fs_folded,
                                        optimizer="log", maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T")


#complex
# Split with continuous symmetric migration, followed by isolation.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "anc_sym_mig", Models_2D.anc_sym_mig, rounds, 5, fs_folded=fs_folded,
                                        optimizer="log", maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T1, T2")


# Split with continuous asymmetric migration, followed by isolation.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "anc_asym_mig", Models_2D.anc_asym_mig, rounds, 6, fs_folded=fs_folded,
                                        optimizer="log", maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, T1, T2")


# Split with no gene flow, followed by period of continuous symmetrical gene flow.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sec_contact_sym_mig", Models_2D.sec_contact_sym_mig, rounds, 5, fs_folded=fs_folded,
                                        optimizer="log", maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T1, T2")


# Split with no gene flow, followed by period of continuous asymmetrical gene flow.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sec_contact_asym_mig", Models_2D.sec_contact_asym_mig, rounds, 6, fs_folded=fs_folded,
                                        optimizer="log", maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, T1, T2")




#Size change
# Split with no migration, then instantaneous size change with no migration.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "no_mig_size", Models_2D.no_mig_size, rounds, 6, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, T1, T2")

# Split with symmetric migration, then instantaneous size change with continuous symmetric migration.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sym_mig_size", Models_2D.sym_mig_size, rounds, 7, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m, T1, T2")


# Split with different migration rates, then instantaneous size change with continuous asymmetric migration.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "asym_mig_size", Models_2D.asym_mig_size, rounds, 8, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2")


