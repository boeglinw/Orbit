#
# plot some orbits , combined top and side view
#
import numpy as np

import get_orbits as go
import get_flux as gf
import get_limiter as gl

import LT.box as B

import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator



plot_em = True
em_filled = True

plot_psirel = True
psirel_filled = False

# set true to plot all orbits
plot_ensemble = False
#plot_ensemble = False

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

#machine = 'NSTX'
machine = 'MAST'

orbit_dir = '/nml_orb_g_MAST_29879_254_g/'
# orbit_dir = '/nml_orb_MAST_p6/'
#orbit_dir = '/nml_orb_MAST_p2_c/'

#orbit_dir = '/nml_orb_p_test/'
#orbit_dir = '/nml_orb_p_lower_position/'
#orbit_dir = '/nml_orb_ppro_5/'
#orbit_dir = '/nml_orb_ppro_6/'

#orbit_dir = '/nml_orb_ppro_7/'
#orbit_dir = '/nml_orb_ppro_7_bundle/'

#orbit_dir = '/nml_orb_ppro_7_high_B/'
#orbit_dir = '/nml_orb_ppro_7_bundle_high_B/'

#orbit_dir = '/nml_orb_p7_upper_position_high_B/'

#orbit_dir = '/nml_orb_p8_upper_position_high_B/'
#orbit_dir = '/nml_orb_p8_lower_position_high_B/'

#orbit_dir = '/nml_orb_p7_lower_position_high_B/'
#orbit_dir = '/nml_orb_p7_bundle_lower_position_high_B/'

#orbit_dir = '/nml_orb_p7_mid_position_high_B/'
#orbit_dir = '/nml_orb_p7_mid_position/'

#orbit_dir = '/nml_orb_p6_lower_position/'
#orbit_dir = '/orb_ppro_5/'
#orbit_dir = '/orb_ppro_8/'
# orbit_dir = '/orb_ppro_6/'

#orbit_dir = './nml_orb_MAST_p0/'
orbit_dir = './' + machine + '_output' + orbit_dir

NC_data_dir = '../Neutron_Camera/'

NC_views = ['view_0.data',\
        'view_1.data',\
        'view_2.data',\
        'view_3.data']


# get the orbits data
NC_tracks = []
for f in NC_views:
    ncf = NC_data_dir + f 
    ncd = B.get_file(ncf)
    nc_r = B.get_data(ncd, 'r')
    nc_x = B.get_data(ncd, 'x')
    nc_y = B.get_data(ncd, 'y')
    nc_z = B.get_data(ncd, 'z')
    NC_tracks.append([nc_r, nc_x, nc_y, nc_z])


#orbit_output = open(orbit_dir + 'orbit_output').readlines()
orbit_output = open(orbit_dir + 'orbit_output').readlines()

# find the EQ file used:
eq_file = 'generic'
for d in orbit_output:
    if (d.find('--> EQ File unit, name') >= 0.):
        eq_file = d.split()[-1:][0]

# flux
print 'reading flux data'
fl = gf.flux(orbit_dir + 'flux.data')
# limiter
print 'reading limiter data'
li = gl.limiter(orbit_dir + 'limiter_drawing.data')
#orbits
print 'reading orbits data'
o = go.orbit(orbit_dir+'orbits.data', fast = True)

# draw side view
f1 = pl.figure(1, figsize= (11,6))
f1.text(0.1, 0.925, eq_file)

# draw 3 regions
li.draw_all()
# draw the rel. flux
# get a nice set of contour lines
# select the first plot
psirel_cont = None
Em_cont = None

axes = li.ax1.get_axes()
f1.sca( axes )
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
# make a nice x-axis
axes.xaxis.set_major_locator(MaxNLocator(4))

# draw a few orbits
if plot_ensemble:
    for i,oo in enumerate(o.orbits):
        o.draw_one(i, ix = 3,  color = 'm')
else:
    n_orb = min(n_orbits_to_plot, len(o.orbits) )
    for i in range(n_orb):
        icol = i % (n_colors-1)
        o.draw_one(i, ix = 3,  color = colors[icol])
# draw NC views
for i,tr in enumerate(NC_tracks):
    pl.plot(tr[0], tr[3], color = colors[i])

# draw  orbits into the top view
axes = li.ax2.get_axes()
f1.sca( axes )
# draw a few orbits
if plot_ensemble:
    for i,oo in enumerate(o.orbits):
        o.draw_one(i, iy = 1,  color = 'm')
else:
    n_orb = min(n_orbits_to_plot, len(o.orbits) )
    for i in range(n_orb):
        icol = i % (n_colors-1)
        o.draw_one(i, iy = 1,  color = colors[icol])

# draw NC views
for i,tr in enumerate(NC_tracks):
    pl.plot(tr[1], tr[2], color = colors[i])

# pl.savefig(orbit_dir+'image.png')
pl.show()

