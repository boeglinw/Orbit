#
# plot some orbits , combined top and side view
#
# this uses the information in the dynamic_input file to assign detector id's and colors them appropriately. It does not use oubit_output as information source
#
# add the mid-plane intersections to get and estimate how they are probed. 
#
import numpy as np
import LT.box as B

import get_flux as gf
import get_limiter as gl

import fileinput as FI

from LT import parameterfile as PF
from LT import pdatafile as DF

import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator

import orbit_view_data as vd
import glob as G

import argparse as AG

dtr = np.pi/180.
#----------------------------------------------------------------------
# decode names and channel numbers
#----------------------------------------------------------------------
def get_names_ch(s):
    vf = s.replace(' ','').split(',')
    name = []
    chan = []
    for v in vf:
        fields = v.split('/')
        v_name = fields[0]
        ch = np.array(fields[1:]).astype(int)
        name.append(v_name)
        chan.append(ch)
    return name, chan


#----------------------------------------------------------------------
# find the first mid-plane crossing of trajectory v
#----------------------------------------------------------------------  
def get_zero_crossing(v):
    is_z = np.where(v.zt<=0.)[0]
    if is_z.shape[0] == 0:
        # no crossing
        return np.nan
    i_n = is_z[0]
    i_p = i_n - 1
    zp = v.zt[i_p]
    rp = v.rt[i_p]
    zn = v.zt[i_n]
    rn = v.rt[i_n]
    m = (zp -zn)/(rp - rn)
    r0 = rn - zn/m
    return r0

#----------------------------------------------------------------------    
# get location of magnetic axis
#----------------------------------------------------------------------  

def get_magnetic_axis(of):
    z_ok = False
    r_ok = False
    for line in FI.input(of):
        if line.find('rmaxis')>=0:
            rmaxis = float(line.split('=')[1])
            r_ok = True
        if line.find('zmaxis')>=0:
            zmaxis = float(line.split('=')[1])
            z_ok = True
        if (r_ok & z_ok):
            FI.close()
            break
    return rmaxis, zmaxis 



#----------------------------------------------------------------------    
# draw orbit in side view
#----------------------------------------------------------------------  
def plot_view_top(PDv, color = 'b', **kwargs):
    first = True
    for v in PDv:
        if first:
            pl.plot(v.xt,v.yt, color = color, **kwargs)
            first = False
        else:
            pl.plot(v.xt,v.yt, color = color)

#----------------------------------------------------------------------  
# draw orbit in top view
#----------------------------------------------------------------------  
def plot_view_side(PDv, color = 'b', **kwargs):
    first = True
    for v in PDv:
        if first:
            pl.plot(v.rt,v.zt, color = color, **kwargs)
            first = False
        else:
            pl.plot(v.rt,v.zt, color = color)


#----------------------------------------------------------------------  
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


all_yellow = pf.get_value('all_yellow',var_type = pf.Bool)

# number of orbits to plot if plot_ensemble is false
n_orbits_to_plot = pf.get_value('n_orbits_to_plot', var_type = int)

# color array
if all_yellow:
    colors = ['y','y','y','y','y','y','y'] # make everything yellow
else:
    colors = ['r','g','b','y','m','k','c','orange']
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

# number of detectors
n_det = pf.get_value('N_det')

# dynamic input file for detector/channels assignment
di_dir = pf.get_value('dynamic_dir')
di_file = pf.get_value('dynamic_file')

dynamic_file = di_dir + '/' + di_file 

dfd = DF.pdfile(dynamic_file)

# get channel/detector assignement
ch,idet = np.array(dfd.get_data_list('ch:detector_id')).T

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

# orbit/Efit data 
orbit_output_dir = orbit_dir

#----------------------------------------------------------------------
# get the detector efficiency from the orbit results
#----------------------------------------------------------------------
# each view can now have several trajectories
PD_views = []
PD_views_f = []
PD_accept = []
PD_channel_number = []
PD_step = []
for i,i_d in enumerate(idet):
    cc = ch[i]
    idd = idet[i]
    name_patt = orbit_output_dir + '/track_{}*.data'.format(i_d)
    PD_view_files = G.glob(name_patt)
    # use the first file in the list to get some parameters used for calculating rates
    PDd = DF.pdfile(PD_view_files[0])
    PD_accept.append(PDd.par.get_value('accept'))
    PD_channel_number.append( PDd.par.get_value('channel_number', var_type = int))
    PD_step.append(PDd.par.get_value('stepsize'))
    # load the trajectories for each view
    PD_v = [ vd.view(f) for f in PD_view_files]
    PD_views_f.append(PD_view_files)
    PD_views.append(PD_v)
    print 'channel : ', cc, ', detecor : ', idd, ' loaded'
PD_accept = np.array(PD_accept)
#----------------------------------------------------------------------
# start drawing
#----------------------------------------------------------------------

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

# axes = li.ax1.get_axes()
# f1.sca( axes )
axes = li.ax1
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

#----------------------------------------------------------------------
# draw orbits
#----------------------------------------------------------------------
# axes = li.ax1.get_axes()
axes = li.ax1
f1.sca( axes )
for i, PDv in enumerate(PD_views):
    icol = ch[i]
    plot_view_side(PDv, color = colors[icol], label = 'Ch {0:d}'.format(ch[i]))

pl.legend(fontsize = 12, loc = 'upper left')
# draw  orbits into the top view
if draw_top_view:
    # axes = li.ax2.get_axes()
    axes = li.ax2
    f1.sca( axes )
    # draw a few orbits
    for i, PDv in enumerate(PD_views):
        icol = ch[i]
        plot_view_top(PDv, color = colors[icol], label = 'Ch {0:d}'.format(ch[i]))
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
    
pl.legend(fontsize = 12, loc = 'upper left')

h_range =(0.6, 1.6)
r0 = []
f2 = pl.figure(2)
for  i,PDv in enumerate( PD_views ):
    icol = ch[i]
    rr = []
    for v in PDv:
        rr.append(get_zero_crossing(v))
    rr = np.array(rr)
    r = rr.mean()
    h = B.histo(rr, range =h_range, bins = 200)
    h.plot(color = colors[icol])
    sig_r = np.sqrt(rr.var())
    r0.append([h, r, sig_r, rr])
# all done
pl.xlabel('R (m)')
pl.ylabel('counts')
pl.title('mid-plane Radius')

# pl.savefig(orbit_dir+'image.png')
pl.show()

def zoom_in(x_lim, y_lim):
    # axes = li.ax1.get_axes()
    axes = li.ax1
    f1.sca( axes )
    pl.xlim(x_lim)
    pl.ylim(y_lim)

def zoom_in2(x_lim, y_lim):
    # axes = li.ax2.get_axes()
    axes = li.ax2
    f1.sca( axes )
    pl.xlim(x_lim)
    pl.ylim(y_lim)


# zoom_in((0.9,1.85),(-.6,.6));zoom_in2((.9,2.),(-.45,.65))

