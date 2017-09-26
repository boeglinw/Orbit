#
# plot the fitted and the simulated source profile
#
import numpy as np

import get_orbits as go
import get_flux as gf
import get_limiter as gl

import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator

import LT.box as B

#----------------------------------------------------------------------
# two gaussians (Input)
#----------------------------------------------------------------------
def Em_2g(r, z):
    # offset gauss in R,z
    g1 = A()*np.exp( -((r-r0())/sig_r())**2)
    g2 = np.exp( -((z-z0())/sig_z())**2)
    return g1*g2    

#----------------------------------------------------------------------
def nearest(xv, x):
    dd = np.abs(x - xv)
    return np.argmin(dd)

plot_em = True
em_filled = True

plot_psirel = True
psirel_filled = False

# set true to plot all orbits
#plot_ensemble = True
plot_ensemble = False

#all_yellow = True
all_yellow = False

# color array
if all_yellow:
    colors = ['y','y','y','y','y','y','y'] # make everything yellow
else:
    colors = ['r','g','b','y','m','k','c']
n_colors = len(colors)
# number of orbits to plot if plot_ensemble is false
n_orbits_to_plot = 8

#orbit_dir = './nml_orb_p_test/'
#orbit_dir = './nml_orb_p_lower_position/'
#orbit_dir = './nml_orb_ppro_5/'
#orbit_dir = './nml_orb_ppro_6/'

#orbit_dir = './nml_orb_ppro_7/'
#orbit_dir = './nml_orb_ppro_7_bundle/'

#orbit_dir = './nml_orb_ppro_7_high_B/'
#orbit_dir = './nml_orb_ppro_7_bundle_high_B/'

#orbit_dir = './nml_orb_p7_upper_position_high_B/'

#orbit_dir = './nml_orb_p8_upper_position_high_B/'
orbit_dir = './nml_orb_p8_lower_position_high_B/'

#orbit_dir = './nml_orb_p7_lower_position_high_B/'
#orbit_dir = './nml_orb_p7_bundle_lower_position_high_B/'

#orbit_dir = './nml_orb_p7_mid_position_high_B/'

#orbit_dir = './nml_orb_p6_lower_position/'
#orbit_dir = './orb_ppro_5/'
#orbit_dir = './orb_ppro_8/'
# orbit_dir = './orb_ppro_6/'

# fit parameters

A  =  B.Parameter(  1.16704568137 )
r0  =  B.Parameter(  1.06498326579 )
sig_r  =  B.Parameter(  0.116676453845 )
z0  =  B.Parameter(  0.0653756066516 )
sig_z  =  B.Parameter(  0.191893739805 )


orbit_output = open(orbit_dir + '/orbit_output').readlines()
# find the EQ file used:
eq_file = 'generic'
for d in orbit_output:
    if (d.find('--> EQ File unit, name') >= 0.):
        eq_file = d.split()[-1:][0]

# flux
print 'reading flux data'
fl = gf.flux(orbit_dir + 'flux.data')
pl.show()

# plot simulated source profile
r_fl = fl.Y[:,0]
z_fl = fl.X[0,:]

z_val = 0.0
r_val = 1.1

r_ind = nearest(r_val, r_fl)
z_ind = nearest(z_val, z_fl)

# the fit function
rr = np.linspace(r_fl.min(), r_fl.max(), 201)
zz = np.linspace(z_fl.min(), z_fl.max(), 201)

Em_r = Em_2g(rr, z_val)
Em_z = Em_2g(r_val, zz)

# plot the line
pl.subplot(1,2,1)
pl.plot(r_fl, fl.Em[:,z_ind], 'ro')
pl.plot(rr, Em_r)

pl.ylim((0., 1.1))
pl.xlim ( (0.5, 1.5) )

pl.xlabel('R (m)')

pl.subplot(1,2,2)
pl.plot(z_fl, fl.Em[r_ind,:], 'ro')
pl.plot(zz, Em_z)

pl.ylim((0., 1.1))
pl.xlim ( (-1., 1.) )

pl.xlabel('z (m)')

# plot contours
pl.figure(2)
pl.subplot(1,2,1)
if plot_em:
    Em_range = fl.Em.max() - fl.Em.min()
    Em_min = fl.Em.min() + Em_range/100.
    v = np.linspace(Em_min, fl.Em.max(), 19)
    if em_filled :
        Em_cont = fl.draw_Em_filled(v)
    else:
        Em_cont = fl.draw_Em(v)

if plot_psirel:
    v = np.linspace(-0.1, fl.Zrel.max(), 23)
    if psirel_filled:
        psirel_cont = fl.draw_psirel_filled(v)
    else:
        psirel_cont = fl.draw_psirel(v)
pl.xlim( (0.8, 1.35) )
pl.ylim( (-0.6, 0.6) )

pl.xlabel('R (m)')
pl.ylabel('z (m)')
pl.title('Simulation')

# now plot fitted function
pl.subplot(1,2,2)
Em_calc = Em_2g(fl.Y, fl.X)
Em_range = Em_calc.max() - Em_calc.min()
Em_min = Em_calc.min() + Em_range/100.
v = np.linspace(Em_min, fl.Em.max(), 19)

pl.contourf(fl.Y, fl.X, Em_calc, v)

if plot_psirel:
    v = np.linspace(-0.1, fl.Zrel.max(), 23)
    if psirel_filled:
        psirel_cont = fl.draw_psirel_filled(v)
    else:
        psirel_cont = fl.draw_psirel(v)

pl.xlim( (0.8, 1.35) )
pl.ylim( (-0.6, 0.6) )

pl.xlabel('R (m)')
pl.title('Fit')

# make a nice x-axis
#axes.xaxis.set_major_locator(MaxNLocator(4))



pl.show()
