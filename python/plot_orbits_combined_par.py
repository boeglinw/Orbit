#
# plot some orbits , combined top and side view
#
import numpy as np

import get_orbits as go
import get_flux as gf
import get_limiter as gl


from LT import parameterfile as PF
from LT import datafile as DF

import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator

import orbit_view_data as vd
import rectangle as RE

import argparse as AG

# contpour plotting controle
cont_scale = 1.0
ncont = 25

draw_top_view = True

# read parameter file
#colormap = pl.get_cmap('CMRmap')
#colormap = pl.get_cmap('gnuplot')
#colormap = pl.get_cmap('gnuplot2')
colormap = pl.get_cmap('gist_heat')
#colormap = pl.get_cmap('jet')
#colormap = pl.get_cmap('hot')

parser = AG.ArgumentParser()
parser.add_argument("control_file", nargs = '?', help="Control file ", default = 'control.data')
args = parser.parse_args()

# open control file
p_file = args.control_file

pf = PF.pfile(p_file)


plot_em = pf.get_value('plot_em',var_type = pf.Bool)
em_filled = pf.get_value('em_filled',var_type = pf.Bool)

plot_psirel = pf.get_value('plot_psirel',var_type = pf.Bool)
psirel_filled = pf.get_value('psirel_filled',var_type = pf.Bool)

# set true to plot all orbits
plot_ensemble = pf.get_value('plot_ensemble',var_type = pf.Bool)
#plot_ensemble = False
if plot_ensemble:
    try:
        ensemble_color = pf.get_value('ensemble_color')
    except:
        ensemble_color = 'm'

#all_yellow = True
all_yellow = pf.get_value('all_yellow',var_type = pf.Bool)

# number of orbits to plot if plot_ensemble is false
n_orbits_to_plot = pf.get_value('n_orbits_to_plot', var_type = int)

# color array
if all_yellow:
    colors = ['y','y','y','y','y','y','y'] # make everything yellow
else:
    colors = ['r','g','b','y','m','k','c', 'purple']
n_colors = len(colors)

machine = pf.get_value('machine')

input_root = pf.get_value('input_root')

input_dir_ext = pf.get_value('input_dir_ext')

input_dir = input_root + machine + input_dir_ext + '/'

input_file = pf.get_value('input_file')

input_name = input_dir + input_file

input_ext = pf.get_value('nml_ext')

input_file_name = input_name + input_ext

output_root = pf.get_value('output_root')

output_dir_ext = pf.get_value('output_dir_ext')


# flux grid
try:
    flux_data_file = pf.get_value('flux_data_file')
except:
    flux_data_file = 'flux.data'
print 'using : ', flux_data_file, ' for flux and Em data' 

# flux limiter
try:
    flux_limiter_file = pf.get_value('flux_limiter_file')
except:
    flux_limiter_file = 'flux_limit.data'
print 'using : ', flux_limiter_file, ' for flux limit data' 


# plot n-flux at mid-plane
try:
    flux_data_file_mp = pf.get_value('flux_data_file_mp')
except:
    flux_data_file_mp = None



orbit_dir = output_root + machine + output_dir_ext + input_file

orbit_output = open(orbit_dir + '/orbit_output').readlines()

# find the EQ file used:
eq_file = 'generic'
for d in orbit_output:
    if (d.find('--> EQ File unit, name') >= 0.):
        eq_file = d.split()[-1:][0]

# flux
print 'reading flux data'
fl = gf.flux(orbit_dir + '/'+flux_data_file)

print 'reading flux limit data'
fll_d = DF.dfile(orbit_dir + '/'+flux_limiter_file)
r_fll = np.array(fll_d.get_data('xlim'))
z_fll = np.array(fll_d.get_data('ylim'))

# limiter
print 'reading limiter data'
li = gl.limiter(orbit_dir + '/limiter_drawing.data')
#orbits
print 'reading orbits data'
o = go.orbit(orbit_dir+'/orbits.data', fast = True)

# draw side view
if draw_top_view:
    f1 = pl.figure(1, figsize= (11,6))
else:
    f1 = pl.figure(1, figsize= (5,8))
f1.text(0.1, 0.925, eq_file)

# draw 3 regions
if draw_top_view:
    li.draw_all()
else:
    li.draw_side_all()
# draw the rel. flux
# get a nice set of contour lines
# select the first plot
psirel_cont = None
Em_cont = None

axes = li.ax1.get_axes()
f1.sca( axes )
# draw the flux limit
pl.plot(r_fll, z_fll, color = 'm', linewidth = 2.)

if plot_em:
    Em_range = fl.Em.max() - fl.Em.min()
    Em_min = fl.Em.min() + Em_range/100.
    v = np.linspace(Em_min, fl.Em.max()*cont_scale, ncont)
    if em_filled :
        Em_cont = fl.draw_Em_filled(v, cmap = colormap)
    else:
        Em_cont = fl.draw_Em(v)
if plot_psirel:
    v = np.linspace(-0.1, fl.Zrel.max(), 23)
    if psirel_filled:
        psirel_cont = fl.draw_psirel_filled(v, cmap = colormap)
    else:
        psirel_cont = fl.draw_psirel(v)
# make a nice x-axis
axes.xaxis.set_major_locator(MaxNLocator(4))

# draw a few orbits
if plot_ensemble:
    for i,oo in enumerate(o.orbits):
        o.draw_one(i, ix = 3,  color = ensemble_color)
else:
    n_orb = min(n_orbits_to_plot, len(o.orbits) )
    for i in range(n_orb):
        icol = i % (n_colors-1)
        o.draw_one(i, ix = 3,  color = colors[icol])

# draw  orbits into the top view

if draw_top_view:
    axes = li.ax2.get_axes()
    f1.sca( axes )
    # draw a few orbits
    if plot_ensemble:
        for i,oo in enumerate(o.orbits):
            o.draw_one(i, iy = 1,  color = ensemble_color)
    else:
        n_orb = min(n_orbits_to_plot, len(o.orbits) )
        for i in range(n_orb):
            icol = i % (n_colors-1)
            o.draw_one(i, iy = 1,  color = colors[icol])

    # get the mid-plane emmissivity
    if plot_em:
        if (flux_data_file_mp != None):
            # load data file
            d = np.load(orbit_dir + '/'+flux_data_file_mp)
            X = d['X']
            Y = d['Y']
            Em_mid = d['EM_mp']
            Em_range = Em_mid.max() - Em_mid.min()
            Em_min = Em_mid.min() + Em_range/100.
            v = np.linspace(Em_min, Em_mid.max()*cont_scale, ncont)
            if em_filled :
                Em_cont = pl.contourf(X,Y,Em_mid, v,cmap=colormap)
            else:
                Em_cont = pl.contour(X,Y,Em_mid, v,cmap=colormap)
    
pl.savefig(orbit_dir+'image.png')
pl.show()

def zoom_in(x_lim, y_lim):
    axes = li.ax1.get_axes()
    f1.sca( axes )
    pl.xlim(x_lim)
    pl.ylim(y_lim)

def zoom_in2(x_lim, y_lim):
    axes = li.ax2.get_axes()
    f1.sca( axes )
    pl.xlim(x_lim)
    pl.ylim(y_lim)

zoom_in((1.6,1.9),(-.2,.2));zoom_in2((0.2,0.6),(1.45,1.9))

class tile:
    def __init__(self, file_name='tilted_tile.data'):
        self.file = file_name
        self.td = DF.dfile(self.file)
        self.r_pos = np.array( self.td.get_data_list('x:y:z') ).T 
        self.r = np.sqrt(self.r_pos[0,:]**2 + self.r_pos[1,:]**2)
        self.rotate(0.)
    
    def rotate(self,phi):
        # rotate all by a common angle (for adjusments)
        Rm = np.matrix([ [np.cos(phi), -np.sin(phi), 0.],\
                         [np.sin(phi),  np.cos(phi), 0.],\
                          [         0.,           0., 1.]])
        self.r_rot = np.squeeze( (Rm*self.r_pos).A )
        
    def plot_mid(self):
        axes = li.ax2.get_axes()
        f1.sca( axes )
        self.xp = np.append(self.r_rot[0,:], self.r_rot[0,0])
        self.yp = np.append(self.r_rot[1,:], self.r_rot[1,0])
        pl.plot(self.xp,self.yp)
        
    def plot_pol(self):
        axes = li.ax1.get_axes()
        f1.sca( axes )
        self.rp = np.append(self.r, self.r[0])
        self.zp = np.append(self.r_rot[2,:], self.r_rot[2,0])
        pl.plot(self.rp,self.zp)

        

