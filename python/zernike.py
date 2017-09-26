# calculate zernike polynomials
# definition used from mathworld.com
#

import numpy as np
from scipy.special import jacobi

# jacobi returns a poly1d
# this could me made more efficient 

class zernike:
    def __init__(self, n, m):
        """
        return a zernike polynomial object of order n,m
        e.g.
        import zerinke as Z

        z42 = Z.zernike(4,2)
        
        value = z42(r, phi)
        """
        # initialize the polynomial
        # if n -  m is odd this is zero
        nm_check = abs(n-m)%2
        if nm_check == 1:
            self.zero = True
            self.value = 0.
            sucess = True
        elif nm_check == 0 :
            self.zero = False
            if m >= 0:
                self.is_even = True
            else:
                self.is_even = False
            # setup radial part
            self.m = np.abs(m)
            self.n = n
            self.norm = (-1)**(0.5*(self.n-self.m))
            alpha = self.m
            beta = 0
            N = (self.n-self.m)/2.
            # get the jacobi polynomial
            self.p_jacobi = jacobi(N, alpha, beta)
        else:
            return None
        # thats's it

    def __call__(self, r, phi):
        if self.zero:
            self.value = 0.
            return 0.
        # evaluate at r, and phi
        # angle part
        if self.is_even:
            self.z_ang = np.cos(self.m* phi)
        else:
            self.z_ang = np.sin(self.m* phi)
        # radial part
        x = 1. - 2.*r**2
        self.z_rad = self.norm * r**self.m * self.p_jacobi(x)
        self.value = self.z_ang*self.z_rad
        return self.value

# that's all
        
            
        
