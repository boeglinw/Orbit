# integrate psi along the path
# use orbit output
# make sure that one obtains the same number if using the same function
# and the same parameters
#
# modification WB March 2015 include acceptance in rate calculation

import numpy as np
import scipy.integrate as integ
import matplotlib.pyplot as pl
#
# use the datafile version that can read the parameters
from LT.pdatafile import pdfile


class view:
    def __init__(self, file, neutron_camera = False):
        # read datafile for the trajectory
        self.d = pdfile(file)
        # get the data
        self.rt = np.array( self.d.get_data('r') )
        self.xt = np.array( self.d.get_data('x') )
        self.yt = np.array( self.d.get_data('y'))
        self.zt = np.array( self.d.get_data('z'))
        self.psi = np.array( self.d.get_data('psi'))
        self.sdlcs = np.array( self.d.get_data('Sdlcs'))
        self.Emt = np.array( self.d.get_data('Em'))
        # one should make sure that psi has a valid value
        # this should be done in orbit not here
        # in_psi_zero = np.where(self.psi == 0.)
        self.psirel = np.array( self.d.get_data('psirel'))
        # self.psirel[in_psi_zero] = 0.
        # indices where psirel is positive
        # self.in_pos_psirel = np.where(self.psirel > 0.)
        self.in_pos_psirel = (self.psirel > 0.)
        self.pos_psirel = self.psirel[self.in_pos_psirel]
        # data from parameter section
        self.acc = self.d.p.get_value('accept', float)
        self.eff_exp = self.d.p.get_value('effic', float)
        self.step = self.d.p.get_value('stepsize', float)
        # total reaction rate
        self.Ntot = self.d.p.get_value('SdV', float)* 4.*np.pi
        # only those parts where psirel is positivie are used
        self.dx = np.diff(self.xt[self.in_pos_psirel])
        self.dy = np.diff(self.yt[self.in_pos_psirel])
        self.dz = np.diff(self.zt[self.in_pos_psirel])
        self.Ema = self.Emt[ self.in_pos_psirel]
        self.dl = np.sqrt(self.dx**2 + self.dy**2 + self.dz**2)
        self.is_NC = neutron_camera
    #end of init

    def get_eff(self, Em, use_all = False, get_rate = False, select = None):
        # use only the those parts of the track where psirel is positive
        # do not use a simpson integration but a simple sum as in orbit
        # self.Sdl = Em(self.pos_psirel[:-1])*self.dl
        # ld = np.arange(len(self.dl)) 
        # simpson integration
        # self.Sdl_int = integ.simps(self.Sdl, x = ld)
        if select == None:
            sel_slice= slice(0,len(self.rt))
        else:
            sel_slice = select
        in_pos_psirel = self.psirel[sel_slice]>0
        pos_psirel = self.psirel[sel_slice][in_pos_psirel]
        if use_all:
            # Em a function of psirel, R and z
            
            r = self.rt[sel_slice][in_pos_psirel]
            z = self.zt[sel_slice][in_pos_psirel]
            self.Sdl = Em(pos_psirel, r, z)*self.step
        else:
            # Em only a function of psirel
            self.Sdl = Em(pos_psirel)*self.step
        # in this case one fits the true rate
        self.Sdl_int = self.Sdl.sum()
        if get_rate:
            # return the rate w/o normalization, include acceptance
            return self.Sdl_int*self.acc
        else:
            # return efficiency
            self.eff = self.Sdl_int*self.acc/self.Ntot
            return self.eff
