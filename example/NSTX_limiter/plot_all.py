import numpy as np
import LT.box as B

import glob as G

dtr = np.pi/180.

plot_side = False

zp = 0.25

file_patt = 'large_R*data'

files = G.glob(file_patt)


files = [\
    'large_R_0.data',\
    'large_R_01.data',\
    'large_R_11.data',\
    'large_R_02.data',\
    'large_R_2.data',\
    'large_R_3.data',\
    'large_R_4.data',\
    'large_R_5.data',\
    'large_R_6.data']

delta_phi = 1.e-3
n_phi = 100.

rl = []
zl = []
phil = []
# large R data
for f in files:
    d = B.get_file(f)
    rl.append(B.get_data(d, 'r'))
    zl.append(B.get_data(d, 'z'))
    phil.append(d.par.get_value('phi'))
rl = np.array(rl)
zl = np.array(zl)
phil.append(phil[0] + 360.)
phil = np.array(phil)
dphil = np.diff(phil)

# small R data
sd = B.get_file('small_R.data')
rs = sd.get_data('r')
zs = sd.get_data('z')
phis = 0.
dphis = 360.
    
if plot_side:
    # plot poloidal projection
    for i, r in enumerate(rl):
        z = zl[i]
        B.plot_line(r,z)
    
    B.pl.xlim((0.,2.))
    B.pl.ylim((-2.,2.))
    
    B.pl.gca().set_aspect('equal')
    B.pl.xlabel = 'R'
    B.pl.ylabel = 'Z'
    B.pl.title = ''
    # inner side
    B.plot_line(rs, zs)

else:
    R_mid = []
    phi_mid = []
    # plot mid-plane at zp
    # find R-values at zp by liner interpolation
    for i, r in enumerate(rl):
        if i< rl.shape[0]:
            j = i
        else:
            j = i - 1 - rl.shape[0]
        z = zl[j]
        phi0 = phil[i]
        phi1 = phil[i] + dphil[i]
        n_phi_val = int(dphil[i]*dtr/delta_phi)
        R_zp0 = np.interp(zp, z, r)
        R_zp1 = np.interp(zp,zl[j], rl[j])
        # constant R part
        phi_val = np.linspace(phi0, phi1, n_phi_val)
        R_val = R_zp0* np.ones_like(phi_val)
        # add new values
        phi_mid += list(phi_val)
        R_mid += list(R_val)
        # constant phi part
        phi_mid += [phi1, phi1]
        R_mid += [R_zp0, R_zp1]
    R_mid = np.array(R_mid)
    phi_mid = np.array(phi_mid)
    # make a polar plot
    fig = B.pl.figure(figsize=(8,8))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)
    ax.plot(phi_mid*dtr, R_mid)
   
B.pl.show()
