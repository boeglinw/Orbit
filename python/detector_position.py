# calculate detector position and directions to be entered into *.inp files
# of orbit

from numpy import *
from math import *

dtr = pi/180


central = input('enter central angle :')*dtr

start = central - input('min angle (below central angle) :')*dtr
stop = central + input('max angle (above central angle) :')*dtr


#step = input('step size :')*dtr
nstep = int(input('number of steps :'))

# angles = arange(start, stop, step)
angles = linspace(start, stop, nstep)

step = diff(angles)[0]

print '----------------------------------------------------------------------'
print 'start = ', start, '(rad)  ', start/dtr, ' (deg)'
print 'stop = ',stop, '(rad)  ', stop/dtr, ' (deg)'
print 'step = ', step, '(rad)  ', step/dtr, ' (deg)'
print '----------------------------------------------------------------------'
print 'for input file'
print start,'\n',stop,'\n', step
print '----------------------------------------------------------------------'
 
print 'angle values : ', angles
