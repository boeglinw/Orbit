#
# plot some orbits  
#
import numpy as np

import get_orbits as go
import get_flux as gf
import get_limiter as gl

import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator

#----------------------------------------------------------------------
# emissivity model
em_par = np.array([1., 11.45])
def Em_mod(f):
    # simple power law profile
    S = em_par[0]*f**em_par[1]
    return S
#----------------------------------------------------------------------

plot_top = False

plot_em = True
em_filled = True

plot_psirel = True
psirel_filled = False


orbit_dir = './orb_pro_52/'
#orbit_dir = './orb_ppro_5/'

# flux
print 'reading flux data'
fl = gf.flux('flux.data')
# limiter
print 'reading limiter data'
li = gl.limiter('limiter_drawing.data')
#orbits
print 'reading orbits data'
o = go.orbit(orbit_dir+'orbits.data', fast = True)

# draw side view
pl.figure(1)
# draw 3 regions
li.draw_side(ireg = 0)
li.draw_side(ireg = 1)
li.draw_side(ireg = 3)
# draw the rel. flux
# get a nice set of contour lines
if plot_em:
    v = np.linspace(0.01, Em_mod(fl.Zrel.max()), 19)
    if em_filled :
        fl.draw_Em_filled(Em_mod, v)
    else:
        fl.draw_Em(Em_mod, v)
if plot_psirel:
    v = np.linspace(0.05, fl.Zrel.max(), 19)
    if psirel_filled:
        fl.draw_psirel_filled(v)
    else:
        fl.draw_psirel(v)
# make a nice x-axis
axes = pl.gca()

axes.xaxis.set_major_locator(MaxNLocator(4))

colors = [\
    '#ff0000',\
        '#ff6600', \
        '#ff9900', \
        '#ffcc00', \
        '#ffff00', \
        '#ccff00', \
        '#99ff00', \
        '#66ff00', \
        '#00ff66', \
        '#00ffcc', \
        '#00ffff', \
        '#00ccff', \
        '#0099ff', \
        '#0066ff', \
        '#0033ff', \
        '#6600ff', \
        '#9900ff', \
        '#cc00ff', \
        '#ff00ff', \
        '#ff0099']

# draw a few orbits
for i in range( o.counter):
    o.draw_one(i, ix = 3, color = colors[i], lw = 1.)
#

if plot_top:
    # draw top view
    pl.figure(2)
# draw 3 regions
    li.draw_top()
# draw a few orbits
    for i in range(o.counter):
        o.draw_one(i, ix = 3, color = colors[i])
#

    



