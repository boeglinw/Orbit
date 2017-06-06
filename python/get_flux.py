# read flux data file an set it up for plotting
from LT.datafile import dfile 
import numpy as np
import matplotlib.pyplot as pl

class flux:
    def __init__(self,file):
        self.d = dfile(file)
        X_set = []
        Y_set = []
        Z_set = []
        Zrel_set = []
        Em_set = []
        X=[]
        Y=[]
        Z=[]
        Zrel = []
        Em = []
        old_index = 0
        for l in self.d.data:
            i = l['i']
            xval = l['z']
            yval = l['r']
            zval = l['psi']
            zvalr = l['psirel']
            Emval = l['Em']
            if i == old_index:
                X_set.append(xval)
                Y_set.append(yval)
                Z_set.append(zval)
                Zrel_set.append(zvalr)
                Em_set.append(Emval)
            elif X_set != []:
                X.append(X_set)
                Y.append(Y_set)
                Z.append(Z_set)
                Zrel.append(Zrel_set)
                Em.append(Em_set)
                X_set = [xval]
                Y_set = [yval]
                Z_set = [zval]
                Zrel_set = [zvalr]
                Em_set = [Emval]
                old_index = i
            else:
                old_index = i
                X_set.append(xval)
                Y_set.append(yval)
                Z_set.append(zval)
                Zrel_set.append(zvalr)
                Em_set.append(Emval)
        # add last data set
        X.append(X_set)
        Y.append(Y_set)
        Z.append(Z_set)
        Zrel.append(Zrel_set)
        Em.append(Em_set)
        self.X = np.array(X)
        self.Y = np.array(Y)
        self.Z = np.array(Z)
        self.Zrel = np.array(Zrel)
        self.Em = np.array(Em)
        self.ncont = 21
        self.v = np.linspace(self.Z.min(), self.Z.max(), self.ncont)
        self.vrel = np.linspace(self.Zrel.min(), self.Zrel.max(), self.ncont)
        # all data have been read

    def calc_Em(self, func):
        self.Em = np.zeros_like(self.Zrel)
        ip = self.Zrel>0.
        self.Em[ip] = func(self.Zrel[ip], self.Y[ip], self.X[ip])

    def set_ncont(self, n):
        self.ncont = n
        self.v = np.linspace(self.Z.min(), self.Z.max(), self.ncont)
        self.vrel = np.linspace(self.Zrel.min(), self.Zrel.max(), self.ncont)
    
    def draw_psi(self, *args, **kwargs):
        # draw the orbits 
        return pl.contour(self.Y, self.X, self.Z, *args, **kwargs)
        #
    def draw_psi_filled(self, *args, **kwargs):
        # draw the orbits 
        return pl.contourf(self.Y, self.X, self.Z, *args, **kwargs)
        #

    def draw_psirel(self, *args, **kwargs):
        # draw the orbits 
        return pl.contour(self.Y, self.X, self.Zrel, *args, **kwargs)
        #
    def draw_psirel_filled(self, *args, **kwargs):
        # draw the orbits 
        return pl.contourf(self.Y, self.X, self.Zrel, *args, **kwargs)
        #
    def draw_Em(self, *args, **kwargs):
        # draw the orbits 
        return pl.contour(self.Y, self.X, self.Em, *args, **kwargs)
        #
    def draw_Em_filled(self, *args, **kwargs):
        # draw the orbits 
        return pl.contourf(self.Y, self.X, self.Em, *args, **kwargs)
        #
        

