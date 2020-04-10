      real*8 function Em(psi, r, z)
      implicit real*8 (a-h, o-z)
      save
c
c calculate emissivity 
c input psi : relative flux
c       r,z : current position in R-Z plane
c
      integer n, pow, gauss, hollow, mod_pow

c parameters to describe profile
      parameter( n = 10)
      real*8 par(n)

c equilibrium information on magnetic axis
      common/mcom/rmaxis,zmaxis


      common /em_par/model, par
      data power/1/, gauss/2/, hollow/3/, gauss2d/4/, mod_pow/5/

      if (psi .lt. 0) then
         Em = 0.
         return
      endif
      
      rr = r - rmaxis
      zr = z - zmaxis
c
c calculate the emissivity for a given values of the flux
c par is an array of parameteters
c
      if (model.eq.power) then
c parameters : 1 amplitude 
c              2 - exponent
         Em = par(1)*psi**par(2)

      elseif (model .eq. gauss) then
c parameters : 1 - amplitude
c              2 - position
c              3 - width
         Em = par(1)*exp( -((rr - par(2))**2+ zr**2)/par(3)**2 )

      elseif (model .eq. mod_pow) then
c parameters: 1 amplitude
c             2 power
c             3 cos ampl
c             4 sin ampl
c             5 cos phase
c             6 sin phase
c polar angle with respect to magnetic axis
         phi = pol_angle(r,z)
         Em = par(1)*(psi**par(2)*(1. + par(3)*cos(phi + par(4)) +
     >        par(5)* sin(phi + par(6)) ) )

      elseif (model .eq. hollow) then
c parameters : 1 - amplitude gauss 1
c              2 - position
c              3 - width, gaus 1
c              4 - amplitude gauss 2
c              5 - width gauss 2
         g1 = par(1)*exp( -((psi-par(2))/par(3))**2)
         g1 = par(4)*exp( -((psi-par(2))/par(5))**2)
         Em = g1 - g2
      elseif (model .eq. gauss2d) then
c parameters : 1 - amplitude gauss 1
c              2 - position in R
c              3 - width in R, 
c              4 - pos. in Z
c              5 - width in Z
c
c FWHM = 2.35 * width
         g1 = par(1)*exp( -0.5*((r-par(2))/par(3))**2)
         g2 = exp( -0.5*((z-par(4))/par(5))**2)
         Em = g1*g2
      endif
      return
      end


c def pol_angle(rx,ry):
c    cphi = rx/np.sqrt(rx**2+ry**2)
c    phi = np.arccos(cphi)
c    phic = np.where(ry>0, phi, 2.*np.pi - phi) # set all values where rx  < 0 to 2pi - phi
c    return phic

      real*8 function pol_angle(r, z)
      implicit real*8 (a-h, o-z)
      data twopi/6.2831853071795862/
      save
      common/mcom/rmaxis,zmaxis

      rr = r - rmaxis
      zr = z - zmaxis

      cphi = min(1.d0, rr/sqrt(rr**2 + zr**2) )
      phi = acos(cphi)
      if (zr .gt. 0.) then
         pol_angle = phi
      else
         pol_angle = twopi - phi
      endif
      return
      end
