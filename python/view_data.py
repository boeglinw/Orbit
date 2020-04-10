# integrate psi along the path
# make sure that one obtains the same number if using the same function
# and the same parameters
import numpy as np
import scipy.integrate as integ
import matplotlib.pyplot as pl
from datafile import dfile


class view:
    def __init__(self, file):
        # read datafile for the trajectory
        self.d = dfile(file)
        # get the data
        self.xt = np.array( self.d.get_data('xt') )
        self.yt = np.array( self.d.get_data('yt'))
        self.psi = np.array( self.d.get_data('psi'))
        self.acc = self.d.get_data('acc')[0]
        self.Ntot = self.d.get_data('Ntot')[0]
        self.eff_exp = self.d.get_data('eff')[0]
        self.dx = np.diff(self.xt)
        self.dy = np.diff(self.yt)
        self.dl = np.sqrt(self.dx**2 + self.dy**2)
    #end of init

    def get_eff(self, Em):
        self.Sdl = Em(self.psi[:-1])*self.dl
        ld = np.arange(len(self.dl)) 
        # simpson integration
        self.Sdl_int = integ.simps(self.Sdl, x = ld)
        self.eff = self.Sdl_int*self.acc/self.Ntot
        return self.eff
