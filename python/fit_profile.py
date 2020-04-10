# giving initial parameters
# fit data with errors
#
# this version is for orbit205 output

import numpy as np
import gen_fit as gf
import matplotlib.pyplot as pl

# special version for orbit data
import orbit_view_data as vd

import zernike as Z

import sys

#----------------------------------------------------------------------
# model function
#----------------------------------------------------------------------

# hollow profile made with a sum of gaussians
# 0: sigma1, 1: sigma2, 2: fract 3: Stot
em_par = np.array([.5, .125, .3, 1.])


def Em_mod(f, par):
    # hollow profile as SOG
    x = (1. - f)
    g1 = np.exp(-(x/par[0])**2)
    g2 = -par[2]*np.exp(-(x/par[1])**2)
    S = par[3]*(g1+g2)
    return S

em_par = np.array([1., 11.45])
def Em_mod(f, par):
    # simple power law profile
    S = par[0]*f**par[1]
    return S



#----------------------------------------------------------------------
# fit functions
#----------------------------------------------------------------------
def Em_pow(x):
    # simple power law
    s = alpha()*x**lam()
    return s

def Em_pol(x):
    # simple polynomial
    s = p0() + p1()*x + p2()*x**2 + p3()*x**3 + p4()*x**4 + p5()*x**5
    return s

#----------------------------------------------------------------------
z20 = Z.zernike(2,0)
z40 = Z.zernike(4,0)
z60 = Z.zernike(6,0)
z80 = Z.zernike(8,0)
z100 = Z.zernike(10,0)

#----------------------------------------------------------------------

def Em_zern(x):
    # zernike polynmials with m = 0
    s = p0() + p1() * z20(x,0.) + p2() * z40(x,0.) + p3() * z60(x,0.) + p4() * z80(x,0.) + p5() * z100(x, 0.)
    return s

#----------------------------------------------------------------------
# the real fitting part is here
#----------------------------------------------------------------------
# set initial values
# total reaction rate
total_rate = 5.e+13
total_rate = 2.e+14
integration_time = 5.e-3

model = 'zer'
use_data_set = 'hollow_crossed'
use_data_set = 'power'

if model == 'pow':
    # simple power law
    lam = gf.Parameter(10)
    alpha = gf.Parameter(1.)
    Em_func = Em_pow
    print 'initial parameters :', lam(), alpha()
    pars = [alpha, lam]
elif model == 'pol':
    # polynomial
    p0 = gf.Parameter(0.)
    p1 = gf.Parameter(0.)
    p2 = gf.Parameter(0.)
    p3 = gf.Parameter(0.)
    p4 = gf.Parameter(0.)
    p5 = gf.Parameter(1.e-12)
    Em_func = Em_pol
    print 'initial parameters :', p0(), p1(), p2(), p3(), p4(), p5()
    pars = [p0, p1, p2, p3, p4, p5]
elif model == 'zer':
    # polynomial
    p0 = gf.Parameter(0.1)
    p1 = gf.Parameter(0.1)
    p2 = gf.Parameter(0.)
    p3 = gf.Parameter(0.)
    p4 = gf.Parameter(0.)
    p5 = gf.Parameter(0.)
    Em_func = Em_zern
    print 'initial parameters :', p0(), p1(), p2(), p3(), p4(), p5()
    pars = [p0, p1, p2, p3, p4, p5]
else:
    print"unknown model", model
    sys.exit()

# define the emissivity model
#----------------------------------------------------------------------


# data files:
data_sets = []
data_dirs = []

view_name = 'power'
view_dir = './orb/'
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

# runs with only hor. positions
view_name = 'hollow_hor'
view_dir = './em_hollow/'
view_files = [ \
    view_dir+'track_h_neg_1.00.data', \
        view_dir+'track_h_neg_1.50.data', \
        view_dir+'track_h_neg_2.00.data', \
        view_dir+'track_h_neg_3.00.data', \
        view_dir+'track_h_pos_0.00.data', \
        view_dir+'track_h_pos_1.00.data', \
        view_dir+'track_h_pos_1.50.data', \
        view_dir+'track_h_pos_2.00.data', \
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
sigma = 1./exp_err

exp_eff = np.zeros_like(exp_counts)
for i,e in enumerate(exp_counts):
    random_factor = 1. + np.random.normal(scale=sigma[i])
    exp_eff[i] = e*random_factor
#
# add an error
# 



# define the fitting function:
def S_f(x): 
    eff = [ views[i].get_eff(Em_func) for i in x]
    eff = np.array(eff)*Ntot
    return eff

stat = gf.fit(S_f, pars ,\
                  exp_eff, \
                  x = xv, \
                  y_err = exp_err, \
                  full_output=1,\
                  ftol = 0.001)
#
# get fit values
fit_eff = stat['fitted values']
r = exp_eff/fit_eff
dr = exp_err/fit_eff
# plot the results:
# first the data with error bars
pl.ioff()
# what the fit produced
pl.figure(1)
# compare the fitted eff. with the exp. one
pl.subplot(2,1,1)
pl.errorbar(xv, exp_eff, yerr = exp_err, marker = 'o', color='r', ls = 'None')
pl.plot(xv, fit_eff, 'bD')

pl.subplot(2,1,2)
pl.errorbar(xv, r, yerr = dr, marker = 'o', color='r', ls = 'None')

pl.figure(2)
# compare the fitted emissivity with the model (input) one
for v in views:
    pl.plot(v.pos_psirel, Em_func(v.pos_psirel), 'bo' )
# what has been put in
for v in views:
    pl.plot(v.pos_psirel, Em_mod(v.pos_psirel, em_par), 'rD') 
pl.show()

print 'Model used : ', model
for i,p_val in enumerate(pars):
    print 'final parameter : ',i, ' = ', p_val()
# that's it
