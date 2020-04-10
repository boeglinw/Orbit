import numpy as np
import matplotlib.pyplot as pl
import h5py as h5


file = "efitOut_26887_200ms.hdf5"

f = h5.File(file)


# get poloidal flux

pf = f['equilibrium/output/profiles2D/poloidalFlux'].value

PF = pf[0]
# y-axis = z
# x-axis = r
PFT = PF.transpose()

# get positions

rV = f['equilibrium/output/profiles2D/rVector'].value
zV = f['equilibrium/output/profiles2D/zVector'].value

# create the grid  for the positions
r,z = np.meshgrid(rV, zV)

# make a contour plot with 30 contours
pl.figure(1, figsize=(8,10) )
pl.contour(r, z, PFT, 30)
ax = pl.gca()
ax.set_aspect('equal')
pl.title(file)
pl.xlabel('R')
pl.ylabel('z')
pl.show()



