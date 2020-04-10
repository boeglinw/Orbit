# giving initial parameters
# fit data with errors
#
# this version is for orbit205 output

import numpy as np
import LT.box as B

import matplotlib.pyplot as pl

# special version for orbit data
import orbit_view_data as vd

# for Zernike pol.
import zernike as Z

# for Chebyshev poly.
import scipy.special as SP

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


#----------------------------------------------------------------------
# model functions
#----------------------------------------------------------------------

# these are the functions used in orbit to generate the data
# use these for tests
em_par_h = np.array([.5, .125, .3, 1.])

def Em_mod_hgauss(f, par = em_par_h):
    # hollow profile as SOG
    x = (1. - f)
    g1 = np.exp(-(x/par[0])**2)
    g2 = -par[2]*np.exp(-(x/par[1])**2)
    S = par[3]*(g1+g2)
    return S

em_par_1 = np.array([1., 11.45])
def Em_mod_1(f, par = em_par_1):
    # simple power law profile
    S = par[0]*f**par[1]
    return S


em_par_4 = [1., 1.0, 0.08, 0.0, 0.15]
def Em_mod_4(r, z, par=em_par_4):
    # offset gauss in R,z
    g1 = par[0]*np.exp( -0.5*((r-par[1])/par[2])**2)
    g2 = np.exp( -0.5*((z-par[3])/par[4])**2)
    return g1*g2    

#----------------------------------------------------------------------
# fit functions
#----------------------------------------------------------------------
# simple power law
#----------------------------------------------------------------------
def Em_pow(x):
    # simple power law
    s = alpha()*x**lam()
    return s

def Em_pow_der(x):
    return np.array([x**lam(), alpha()*np.log(x)*x**lam()])
#    return np.array([alpha()*np.log(x)*x**lam()])
#----------------------------------------------------------------------
# polynomial
#----------------------------------------------------------------------

def Em_pol(x):
    # simple polynomial
    s = p0() + p1()*x + p2()*x**2 + p3()*x**3 + p4()*x**4 + p5()*x**5
    return s

def Em_pol_der(x):
    return np.array([np.ones_like(x), x, x**2, x**3, x**4, x**5])
#----------------------------------------------------------------------
# polynomial squared
#----------------------------------------------------------------------

def Em_pol_sq(x):
    # simple polynomial
    s = (p0() + p1()*x + p2()*x**2 + p3()*x**3 )**2
    return s

def Em_pol_sq_der(x):
    s = (p0() + p1()*x + p2()*x**2 + p3()*x**3)
    return np.array([np.ones_like(x), x, x**2, x**3])*2.*s

#----------------------------------------------------------------------
z20 = Z.zernike(2,0)
z40 = Z.zernike(4,0)
z60 = Z.zernike(6,0)
z80 = Z.zernike(8,0)
z100 = Z.zernike(10,0)

#----------------------------------------------------------------------
# Zernike polynomials
#----------------------------------------------------------------------

def Em_zern(x):
    # zernike polynmials with m = 0
    s = p0() + p1() * z20(x,0.) + p2() * z40(x,0.) + p3() * z60(x,0.) + p4() * z80(x,0.) + p5() * z100(x, 0.)
    return s

def Em_zern_der(x):
    return np.array([np.ones_like(x),z20(x,0.),z40(x,0.),z60(x,0.),z80(x,0.),z100(x, 0.)])

#----------------------------------------------------------------------
# Chenyshev pol. squared
#----------------------------------------------------------------------
def cheb_sum(x):
    s = p0()*SP.eval_chebyt(0,x) +\
        p1()*SP.eval_chebyt(1,x) +\
        p2()*SP.eval_chebyt(2,x) +\
        p3()*SP.eval_chebyt(3,x) +\
        p4()*SP.eval_chebyt(4,x) 
    return s    
def Em_cheb(x):
    # Chebyshev polynomial
    s = ( cheb_sum(x) ) ** 2
    return s

def Em_cheb_der(x):
    s = np.sign(cheb_sum(x) )
    return 2.*s*np.array([\
        SP.eval_chebyt(0,x) ,\
        SP.eval_chebyt(1,x) ,\
        SP.eval_chebyt(2,x) ,\
        SP.eval_chebyt(3,x) ,\
        SP.eval_chebyt(4,x)  ])

#----------------------------------------------------------------------
# power law with modulation (requires psi, r, z)
#----------------------------------------------------------------------
# from orbit output: rmaxis = 1.071290730999
#                    zmaxis = -1.46114150999999997E-004
r_maxis = 1.071290730999
z_maxis = -1.46114150999999997E-004

def Em_pow_mod(psi, r, z ):
    # r_maxis : location of magnetic axis
    rho_x = r - r_maxis
    rho_z = z - z_maxis
    phi = pol_angle(r,z)  # polar angle
    s = alpha()*psi**lam()*(1.+psi**lam1() * (A1()*np.sin(phi) + B1()*np.cos(phi))**2 )
    return s

def Em_pow_mod_der(psi, r, z ):
    return np.ones_like(r)*0.1

                                          
#----------------------------------------------------------------------
# two gaussians (Input)
#----------------------------------------------------------------------
def Em_2g(psi, r, z):
    # offset gauss in R,z
    g1 = A()*np.exp( -((r-posr())/sigr())**2)
    g2 = np.exp( -((z-posz())/sigz())**2)
    return g1*g2    

def Em_2g_der(psi, r, z):
    # offset gauss in R,z
    return np.ones_like(r) * 0.1

#----------------------------------------------------------------------
# the real fitting part is here
#----------------------------------------------------------------------
# set initial values
# total reaction rate
total_rate = 2.e+14
total_rate = 5.e+13
integration_time = 1e-3

# Default setting is that Em depends only on psirel
use_all_variables = False

model = 'zer'
model = 'cheb'
model = 'pol'
model = 'pow'
model = 'two_gauss'
model = 'pow_mod'

model = 'pow'

use_data_set = 'p7'

if model == 'pow':
    # simple power law
    lam = B.Parameter(11., 'lam')
    alpha = B.Parameter(1., 'alpha')
    Em_func = Em_pow
    Em_func_der = Em_pow_der
    print 'initial parameters :', alpha,'\n', lam,'\n'
    current_fit_par = [alpha, lam]
elif model == 'pol':
    # polynomial
    p0 = B.Parameter(.0)
    p1 = B.Parameter(0.1)
    p2 = B.Parameter(0.1)
    p3 = B.Parameter(1)
#    p4 = B.Parameter(.1)
#    p5 = B.Parameter(.1)
#    Em_func = Em_pol
#    Em_func_der = Em_pol_der
    Em_func = Em_pol_sq
    Em_func_der = Em_pol_sq_der
#    print 'initial parameters :', p0(), p1(), p2(), p3(), p4(), p5()
#    current_fit_par = [p0, p1, p2, p3, p4, p5]
    print 'initial parameters :', p0(), p1(), p2(), p3()
    current_fit_par = [p0, p1, p2, p3]
elif model == 'zer':
    # polynomial
    p0 = B.Parameter(0.1)
    p1 = B.Parameter(0.1)
    p2 = B.Parameter(0.)
    p3 = B.Parameter(0.)
    p4 = B.Parameter(0.)
    p5 = B.Parameter(0.)
    Em_func = Em_zern
    Em_func_der = Em_zern_der
    print 'initial parameters :', p0(), p1(), p2(), p3(), p4(), p5()
    current_fit_par = [p0, p1, p2, p3, p4, p5]
elif model == 'cheb':
    # polynomial
    p0 = B.Parameter(0.1)
    p1 = B.Parameter(0.1)
    p2 = B.Parameter(0.)
    p3 = B.Parameter(0.)
    p4 = B.Parameter(0.)
#    p5 = B.Parameter(0.)
    Em_func = Em_cheb
    Em_func_der = Em_cheb_der
    print 'initial parameters :', p0(), p1(), p2(), p3(), p4()# , p4(), p5()
    current_fit_par = [p0, p1, p2, p3, p4]
elif model == 'pow_mod':
    # modulated power law
    alpha = B.Parameter(.8, 'alpha')
    lam   = B.Parameter(15., 'lam')
    lam1   = B.Parameter(0., 'lam1')
    A1    = B.Parameter(0.2, 'A1')
    B1    = B.Parameter(0.2, 'B1')
    A2    = B.Parameter(1.001)
    B2    = B.Parameter(1.001)
    use_all_variables = True
    Em_func = Em_pow_mod
    Em_func_der = Em_pow_mod_der
    print 'initial parameters :', alpha(), lam(), A1(), B1(), A2(), B2()
    current_fit_par = [alpha, lam, A1, B1]
elif model == 'two_gauss':
    # modulated power law
    A = B.Parameter(1., 'A')
    posr   = B.Parameter(1.,'r0')
    sigr    = B.Parameter(0.1,'sig_r')
    posz   = B.Parameter(0.1,'z0')
    sigz    = B.Parameter(0.1,'sig_z')
    use_all_variables = True
    Em_func = Em_2g
    Em_func_der = Em_2g_der
    print 'initial parameters :', A,'\n', posr,'\n', sigr,'\n', posz,'\n', sigz, '\n'
    current_fit_par = [A, posr, sigr, posz, sigz]
else:
    print"unknown model", model
    sys.exit()

# define the emissivity model
#----------------------------------------------------------------------


# data files:
data_sets = []
data_dirs = []

view_name = 'power'
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


data_sets = {view_name: view_files}

# data sets for p7_upper_position and p7_lower_position 
view_name = 'p7'
view_dir = './nml_orb_p7_upper_position_high_B/'
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
view_dir = './nml_orb_p7_lower_position_high_B/'
view_files += [ \
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

# define the fitting function:
def S_f(x): 
    eff = [ views[i].get_eff(Em_func, use_all = use_all_variables) for i in x]
    eff = np.array(eff)*Ntot
    return eff

orbit_fit = B.genfit(S_f, current_fit_par ,\
              y = exp_eff, \
              x = xv, \
              y_err = exp_err, \
              full_output=1,\
                nplot = 0, \
              ftol = 0.001, maxfev = 2000
                  )
# fitting using analytically calculated derivatives does not work well
#                  deriv = Em_func_der, \
#
# get fit values
stat = orbit_fit.stat
fit_eff = stat['fitted values']
# get covariance matrix
if orbit_fit.covar == None:
    print "fit did not converge !' "
    for k in stat.keys():
        print '------------------------------------------------------'
        print k, ' = ', stat[k]
    print '------------------------------------------------------'
    sys.exit()
#
# take chi square into account for error estimate
mcov = orbit_fit.covar *  orbit_fit.chi2_red

# ratio between exp. and fit
r = exp_eff/fit_eff
dr = exp_err/fit_eff

# plot the results:
# first the data with error bars
pl.ioff()
# what the fit produced
fig1 = pl.figure(1, figsize = figure_size_2)
# compare the fitted eff. with the exp. one
pl.subplot(2,1,1)
pl.errorbar(xv, exp_eff, yerr = exp_err, marker = 'o', color='r', ls = 'None')
pl.plot(xv, fit_eff, 'bD')
pl.ylabel('counts')
ymin,ymax = pl.ylim()
if ymin < 0.:
    pl.ylim( (0., ymax) )
# done
pl.subplot(2,1,2)
# compare the fitted emissivity with the model (input) one
psirel = np.linspace(0., 1., 101)
# calculate errors for the fittet emissivity at each psi value
# using the uncertainty of the fit
if use_all_variables:
    for i,vv in enumerate(views):
        ic = i%7
        Em_fit_view = Em_func(vv.pos_psirel, vv.rt[vv.in_pos_psirel], vv.zt[vv.in_pos_psirel])
        Em_diff_view = (Em_fit_view - vv.Ema)/vv.Ema* 100.
        pl.plot(vv.pos_psirel, Em_fit_view, color = colors[ic] )
        pl.plot(vv.pos_psirel, vv.Ema, color = colors[ic], ls = '--' )
else:
    dEm_dpar = Em_func_der(psirel)
    der = np.matrix( dEm_dpar )
    Em_fit = Em_func(psirel)
    ss=(der.transpose()*mcov)*der
    # these are the errors for each view
    de_err = np.sqrt( np.diag(ss) )
    # fitted function
    Em_high = Em_fit + de_err
    Em_low =  Em_fit - de_err
    pl.fill_between(psirel, Em_high, Em_low, color = 'b')
    pl.plot(psirel, Em_fit, color = 'c', lw = 2.)
# use the model from Orbit
#for vv in views:
#    pl.plot(vv.pos_psirel, vv.Ema)
#
pl.xlabel(r'$\psi_{rel}$')
pl.ylabel(r'$S(\psi_{rel})$')
# done
pl.subplots_adjust(left=0.20, bottom = 0.12)

fig2 = pl.figure(2, figsize = figure_size)
pl.errorbar(xv, r, yerr = dr, marker = 'o', color='r', ls = 'None')
ymin,ymax = pl.ylim()
pl.xlabel('view number')
pl.ylabel('ratio exp/fit')
if ymin < 0.8:
    pl.ylim( (0.8, ymax) )
ymin,ymax = pl.ylim()
if ymax > 1.2:
    pl.ylim( (ymin, 1.2) )

pl.subplots_adjust(left=0.15, bottom = 0.15)

print 'reduced chi square: ', orbit_fit.chi2_red
print 'Model used : ', model
for i,p_val in enumerate(current_fit_par):
    print 'final parameter : ',i, ' = ', p_val
# that's it
