from LT.datafile import dfile 
import numpy as np
import matplotlib.pyplot as pl

class orbit:
    def __init__(self,file, fast = True):
        # using the fast option does not use the datafile module
        # it is much faster that way but less general
        if fast:
            self.data = open(file).readlines()[1:]
        else:
            self.d = dfile(file)
        self.orbits = []
        new_orbit = []
        counter = 0
        if fast :
            for l in self.data:
                f = map(float, l.split())
                x = f[1]*np.cos( f[2] )
                y = f[1]*np.sin( f[2] )
                # this creates an x, y, z, r, phi, z, ... array
                # f[1:] to skip the step number
                r = [x, y, f[3] ] + f[1:]
                data = np.array(r)
                if f[0] == 1:
                    counter += 1
                    print "add orbit ", counter
                    if new_orbit != []:
                        self.orbits.append(np.vstack(new_orbit))
                        new_orbit = []
                new_orbit.append(data)
        else:
            for l in self.d.data:
                dt = ([\
                        l['r'], l['phi'],\
                            l['z'], l['vr'], l['vphi'], l['vz'],\
                            l['br'], l['bphi'], l['bz']\
                            ])
                x = l['r']*np.cos( l['phi'] )
                y = l['r']*np.sin( l['phi'] )
                # this creates an x, y, z, r, phi, z, ... array
                r = [x, y, l['z'] ]
                data = np.array(r + dt)
                if l['step'] == 1:
                    counter += 1
                    print "add orbit ", counter
                    if new_orbit != []:
                        self.orbits.append(np.vstack(new_orbit))
                        new_orbit = []
                new_orbit.append(data)
        # add the last orbit
        counter += 1
        print "add last orbit "
        self.orbits.append(np.vstack(new_orbit))
        # all data have been read
        self.counter = counter-1
        print "total of : ", self.counter, " orbits loaded !"
    
    def draw(self, ix = 0, iy = 2, *args, **kwargs):
        # draw the orbits 
        for o in self.orbits:
            pl.plot(o[:,ix], o[:,iy], *args, **kwargs)
        #
    def is_ok(self, i):
        # check the orbit index
        if i < 0:
            print 'invalid orbit index : ', i
            return False
        if i > (len(self.orbits)-1):
            print 'orbit index too large :',i,  ' max = ', len(self.orbits)-1
            return False
        return True

    def draw_one(self, iorb, ix = 0, iy = 2, *args, **kwargs):
        # draw the orbits 
        if self.is_ok(iorb):
            pl.plot(self.orbits[iorb][:,ix], self.orbits[iorb][:,iy], *args, **kwargs)
        #

        
        

