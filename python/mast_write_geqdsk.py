import numpy as np
import matplotlib.pyplot as pl
import itertools

import h5py as h5

# defin grouping iterator
def group(lst, n):
      for i in range(0, len(lst), n):
            val = lst[i:i+n]
            yield tuple(val)

def group3(lst, n, defaultvalue = 0):
        """group([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 3) => iterator
    
        Group an iterable into an n-tuples iterable. Incomplete tuples
        are padded with defaultvalue
    
        >>> list(group(range(10), 3))
        [(0, 1, 2), (3, 4, 5), (6, 7, 8), (9, 0, 0)]
        """
        # array of slices starting at range(n) and stepping over n elements
        # first element: [0, 3, 6, 9]
        # second       : [1, 4, 7]
        # third        : [2, 5, 8]
        slices_array = [itertools.islice(lst, i, None, n) for i in range(n)]
        # zip the arrays together fill missing elements with fill value
        # *slices_array passes the arrays as an arbitrary number of argument
        return itertools.izip_longest( *slices_array, fillvalue = defaultvalue)


def print_group5(x, out):
    for g in group3(x, 5):
        out.write(std_line_2020(g[0], g[1], g[2], g[3], g[4]))


# Formatting according to fortran

def std_line_2000(s0,i1,i2,i3):
    return "{0:" "<48s}{1:4d}{2:4d}{3:4d}\n".format(s0, int(i1), int(i2), int(i3))

def std_line_2020(x0,x1,x2,x3,x4):
    return "{0:16.9e}{1:16.9e}{2:16.9e}{3:16.9e}{4:16.9e}\n".format(x0,x1,x2,x3,x4)

def std_line_2022(x0,x1):
    return "{0:5d}{1:5d}\n".format(x0,x1)



#file = "efitOut_26887_200ms.hdf5"
file = "efitOut_27148.hdf5"

f = h5.File(file)


# header information
# read (neqdsk,2000) (case(i),i=1,6),idum,nw,nh

pulse_number = int(f['/equilibrium/header/pulseNumber'].value[0])
times = f['/equilibrium/header/times'].value[0]*1000 #in ms

# hor grid points
nw = f['/equilibrium/input/regularGrid/nr'].value[0]
# ver grod points
nz = f['/equilibrium/input/regularGrid/nz'].value[0]

# read (neqdsk,2020) rdim,zdim,rcentr,rleft,zmid
rmin = f['/equilibrium/input/regularGrid/rMin'].value[0]
rmax = f['/equilibrium/input/regularGrid/rMax'].value[0]

rdim = rmax - rmin
rcentr = 0.5*(rmax + rmin)
rleft = rmin


zmin = f['/equilibrium/input/regularGrid/zMin'].value[0]
zmax = f['/equilibrium/input/regularGrid/zMax'].value[0]
zmid = 0.5*(zmin+zmax)
zdim = zmax - zmin

simag = f['/equilibrium/output/globalParameters/psiAxis'].value[0]
sibry = f['/equilibrium/output/globalParameters/psiBoundary'].value[0]


# magnetic axis
rmaxis,zmaxis=f['/equilibrium/output/globalParameters/magneticAxis'].value[0]
# magnetic field there
bvacRmag = f['/equilibrium/output/globalParameters/bvacRmag'].value[0]

# geometric axis
rgeom,zgeom=f['/equilibrium/output/geometry/geometricAxis'].value[0]
# magnetic field there
bvacRgeom = f['/equilibrium/output/globalParameters/bvacRgeom'].value[0]



#plasma current
current = f['/equilibrium/output/globalParameters/plasmaCurrent'].value[0]

fpol = f['/equilibrium/output/fluxFunctionProfiles/fDia'].value[0]
ffprime = f['/equilibrium/output/fluxFunctionProfiles/ffprime'].value[0]

pres = f['/equilibrium/output/fluxFunctionProfiles/staticPressure'].value[0]
pprime = f['/equilibrium/output/fluxFunctionProfiles/staticPprime'].value[0]

qpsi = f['/equilibrium/output/fluxFunctionProfiles/q'].value[0]

# get poloidal flux

psirz = f['equilibrium/output/profiles2D/poloidalFlux'].value[0]
psirzt = psirz.transpose()

# get positions

rV = f['equilibrium/output/profiles2D/rVector'].value
zV = f['equilibrium/output/profiles2D/zVector'].value

# boundary
bndry=f['/equilibrium/output/geometry/boundaryCoordinates'].value[0]
# r of boundary
r_bndry = bndry[:,0]
# z of boundary
z_bndry = bndry[:,1]
nbbbs = len(r_bndry)

# limiter values
rlim = f['/equilibrium/input/limiter/rValues'].value
zlim = f['/equilibrium/input/limiter/zValues'].value
limitr = len(zlim)

idum = 0
xdum = 0.


o = open("gMAST_shot_{0:5d}.dat".format(pulse_number), 'w')

# print EQDSK values
#print "MAST shot: ", pulse_number, ", time = ", times, "ms", idum, nw, nz
title = "MAST shot: {0:5d}, time = {1:7.1f}".format(pulse_number, times)
o.write(std_line_2000(title, idum, nw, nz))
#print rdim, zdim, rgeom, rleft, zmid
o.write(std_line_2020(rdim, zdim, rgeom, rleft, zmid))

#print rmaxis, zmaxis, simag, sibry, bvacRgeom
o.write(std_line_2020(rmaxis, zmaxis, simag, sibry, bvacRgeom))
#print current, simag, xdum, rmaxis, xdum
o.write(std_line_2020(current, simag, xdum, rmaxis, xdum))
#print zmaxis, xdum, sibry, xdum, xdum
o.write(std_line_2020(zmaxis, xdum, sibry, xdum, xdum))
#print fpol
print_group5(fpol, o)
#print pres
print_group5(pres, o)
#print ffprime
print_group5(ffprime, o)
#print pprime
print_group5(pprime, o)
#print psirz
# might need to take the transpose()
print_group5(psirzt.flatten(), o)
#print qpsi
print_group5(qpsi, o)
#print nbbbs, limitr
o.write( std_line_2022(nbbbs, limitr) )
# prepare the boundaries
rz_bndry = ( np.vstack((r_bndry,z_bndry)) ).transpose()
print_group5(rz_bndry.flatten(), o)

rz_lim =  ( np.vstack((rlim,zlim)) ).transpose()
print_group5(rz_lim.flatten(), o)
# all done
    
o.close()

