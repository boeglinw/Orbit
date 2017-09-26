# fir random gaussians to orbit data
#
# this version is for orbit205 output

import numpy as np
import LT.box as B

import matplotlib.pyplot as pl

# special version for orbit data
import orbit_view_data as vd

import pdb

import sys

figure_size = (8,6)
figure_size_2 = (6,8)

colors = ['r','g','b','y','c','m','k']

#----------------------------------------------------------------------
# clear all figures
def cl_all(n):
    for i in range(1,n+1):
        pl.figure(i);pl.clf();
    # that is all
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# get a polar angle for a 2d vector
def pol_angle(rx,ry):
    cphi = rx/np.sqrt(rx**2+ry**2)
    phi = np.arccos(cphi)
    phic = np.where(ry>0, phi, 2.*np.pi - phi) # set all values where rx  < 0 to 2pi - phi
    return phic
#----------------------------------------------------------------------

def print_init_par(par):
    print "Initial parameters :"
    for i,p in enumerate(par):
        print "parameter : ", i, " = ", p


#----------------------------------------------------------------------
# 2d-gaussian function class
#----------------------------------------------------------------------
class EM:
    def __init__(self,r0=0., sigr=1., z0=0., sigz=1):
        # offset gauss in R,z
        self.r0 = r0
        self.sigr = sigr
        self.z0 = z0
        self.sigz = sigz

    def set_par(self, r0=0., sigr=1., z0=0., sigz=1):
        # offset gauss in R,z
        self.r0 = r0
        self.sigr = sigr
        self.z0 = z0
        self.sigz = sigz

    def gauss(self,psirel, r,z):   
        gr = np.exp( -0.5*((r-self.r0)/self.sigr)**2)
        gz = np.exp( -0.5*((z-self.z0)/self.sigz)**2)
        return gr*gz

    def __call__(self,x):
        # eff = [ views[i].get_eff(self.gauss, use_all = use_all_variables) for i in x]
        eff = views[int(x)].get_eff(self.gauss, use_all = use_all_variables)
        return eff

#----------------------------------------------------------------------
# fit function parameters
#----------------------------------------------------------------------
# total number of gaussians
Ngauss = 4
# range of sigma in R
sigr_min = 0.02
sigr_max = 0.2
sigr_range = sigr_max - sigr_min

# range of sigma in z
sigz_min = 0.05
sigz_max = 0.2
sigz_range = sigz_max - sigz_min

# from orbit output: magnetic axis
#                    rmaxis = 1.071290730999
#                    zmaxis = -1.46114150999999997E-004
r_maxis = 1.071290730999
z_maxis = -1.46114150999999997E-004

# position ranges
r_min = r_maxis - 0.3
r_max = r_maxis + 0.3
r_range = r_max - r_min

z_min = z_maxis - 0.3
z_max = z_maxis + 0.3
z_range = z_max - z_min

# number of trials
N_fits = 5000


#----------------------------------------------------------------------
# the real fitting part is here
#----------------------------------------------------------------------
# set initial values
# total reaction rate
total_rate = 5.e+13
total_rate = 2.e+14
integration_time = 5e-3

# Default setting is that Em depends only on psirel
use_all_variables = True

use_data_set = '2gauss'



#----------------------------------------------------------------------
# setup the data files
#----------------------------------------------------------------------
# data files:
data_sets = []
data_dirs = []

# 2 gauss model
view_name = '2gauss'
# set1
view_dir = './nml_orb_ppro_5/'
view_files = [ \
    view_dir+'track_1011.data', \
    view_dir+'track_2011.data', \
    view_dir+'track_3011.data', \
    view_dir+'track_4011.data', \
    view_dir+'track_5011.data', \
    view_dir+'track_6011.data', \
    view_dir+'track_7011.data', \
    view_dir+'track_8011.data' \
    ]
# add set 2
view_dir = './nml_orb_p_lower_position/'
view_files += [ \
    view_dir+'track_1011.data', \
    view_dir+'track_2011.data', \
    view_dir+'track_3011.data', \
    view_dir+'track_4011.data', \
    view_dir+'track_5011.data', \
    view_dir+'track_6011.data', \
    view_dir+'track_7011.data', \
#   view_dir+'track_8011.data' \
]

data_sets = {view_name: view_files}

# runs with hor and vert. positions
view_name = 'hollow_crossed'
view_dir = './em_hollow_crossed/'
view_files = [ \
    view_dir+'track_h_neg_1.00.data', \
        view_dir+'track_h_neg_2.00.data', \
        view_dir+'track_h_pos_0.00.data', \
        view_dir+'track_h_pos_1.00.data', \
        view_dir+'track_v_neg_1.00.data', \
        view_dir+'track_v_pos_0.00.data', \
        view_dir+'track_v_neg_2.00.data', \
        view_dir+'track_v_pos_1.00.data', \
        ]
data_sets[view_name]= view_files

# initialize data and setup exp. data file
view_files = data_sets[use_data_set]
Ntot = total_rate * integration_time

views = [ vd.view(f) for f in view_files]
exp_dat = [view.eff_exp for view in views ]
exp_counts = np.array(exp_dat)*Ntot
xv = np.arange(len(views))


exp_err = np.sqrt(exp_counts)
exp_err.clip(1., out=exp_err)*1.  # make sure the min.error is 1 not 0

sigma = 1./exp_err

exp_eff = np.zeros_like(exp_counts)
for i,e in enumerate(exp_counts):
    random_factor = 1. + np.random.normal(scale=sigma[i])
    exp_eff[i] = max(e * random_factor, 0.)
#
# add an error
#
penn_factor = 1.e30

#----------------------------------------------------------------------
# do the fit
#----------------------------------------------------------------------
# list of good fits
good_fit = False
fit_results = []
#
k = 0
cc = 0
while k <= N_fits:
    cc += 1
    if cc > 10000:
        print 'exhausted 1000 tries !'
        break
    # define gaussians
    func = []
    print 'fit : ', k
    for ng in range( Ngauss):
        # select widths
        sigr = np.random.uniform()*sigr_range + sigr_min
        sigz = np.random.uniform()*sigz_range + sigz_min
        #select positions
        rg = r_min + np.random.uniform()*r_range
        zg = z_min + np.random.uniform()*z_range
        func.append( EM(rg, sigr, zg, sigz) )
    # now do the fit
    try:
        orbit_fit = B.gen_linfit(func, \
                                 y = exp_eff, \
                                 x = xv, \
                                 yerr = exp_err )
        good_fit = True
    except:
        print 'problem with fit '
        good_fit = False
    # get the chi-square
    if good_fit and (orbit_fit.chi_red < 500.):
        k += 1
        fit_results.append([orbit_fit.chi_red, orbit_fit.par, func])
# finished with all possible fits
fit_results = np.array(fit_results)

def get_em(r,z, par, func):
    sum = 0.
    psi = 0.
    for i,f in enumerate(func):
        sum += par[i] * f.gauss(psi, r, z)
    return sum

def plot_em_r(rmin, rmax, z, chi_red_min):
    rr = np.linspace( rmin, rmax, 100)
    for f in fit_results:
        if (f[0] <= chi_red_min):
            em = get_em(rr, z, f[1], f[2])
            plot(rr, em, 'b')

def plot_em_z(zmin, zmax, r, chi_red_min):
    zz = np.linspace( zmin, zmax, 100)
    for f in fit_results:
        if (f[0] <= chi_red_min):
            em = get_em(r, zz, f[1], f[2])
            plot(zz, em, 'r')
    


