import numpy as np

class EM:
    def __init__(self,r0=0., sigr=1., z0=0., sigz=1):
        # offset gauss in R,z
        self.r0 = r0
        self.sigr = sigr
        self.z0 = z0
        self.sigz = sigz

    def set_par(self, r0=0., sigr=1., z0=0., sigz=1):
        # offset gauss in R,z
        self.r0 = r0
        self.sigr = sigr
        self.z0 = z0
        self.sigz = sigz

    def __call__(self,r,z):
        gr = np.exp( -0.5*((r-self.r0)/self.sigr)**2)
        gz = np.exp( -0.5*((z-self.z0)/self.sigz)**2)
        return gr*gz
