

       subroutine intlim (hitlim)
c************************************************************
c This routine checks to see whether a particle has struck the
c limiter on its most recent step in its orbit. Return the 
c logical variable HITLIM as TRUE if particle has hit the 
c limiter, FALSE otherwise.  The particle is considered to have
c hit the limiter if it is on the line forming the outline of
c the limiter, or is on the side of that line away from where
c the plasma is supposed to reside.
c************************************************************

      include 'dcalc.h'

      save
      

c      r assumed to be the R, phi, and z coordinates, of particle
c		plus total distance from the origin

      parameter (ntorrmx = 30)
      parameter (nlplmx = 80)
      parameter (nlprmx = 80)

c     NTORRMX is the maximum number of toroidal regions (ie intervals in
c       phi, the toroidal angle).  Used to dimension arrays.
c     NLPLMX is the maximum number of points for defining the left part
c       (ie inboard part) of the limiter
c     NLPRMX maximum number of points for right (outboard) part of limiter


      common /lim/ ntorreg, phistt(ntorrmx+1), nlimplf(ntorrmx),
     .  nlimprt(ntorrmx), rlimplf(nlplmx, ntorrmx), 
     .  zlimplf(nlplmx, ntorrmx), rlimprt(nlprmx, ntorrmx),
     .  zlimprt(nlprmx, ntorrmx)
      integer ntorreg, nlimplf, nlimprt
      real*8 phistt, rlimplf, zlimplf, rlimprt, zlimptrt

c     NTORREG number of toroidal regions (ie intervals in phi, the toroidal
c       angle).  Must be .LE. NTORRMX
c     PHISTT toroidal angle at which each toroidal region starts (ending
c       angle is the same as the starting angle of the next region, and
c       there is one extra element at the end of this array, which is 
c       set, by definition, to 2*pi)
c     NLIMPLF, NLIMPRT number of points describing, respectively, the left
c       and right halves of limiter for this particular toroidal region.
c     RLIMPLF, ZLIMPLF R & z coordinates of points which define the limiter
c       on the left side (small R side) of the plasma for each toroidal 
c       region
c     RLIMPRT, ZLIMPRT same as above, but for right (large R) side of plasma

      common /consts/ pm, echrg, emass, pi
      parameter (tupi = 6.283185307180)
      
       logical hitlim
       integer srchz
	   
	   
c set HITLIM to .TRUE. if particle has hit limiter.
       hitlim = .false.


c First, compute phi (toroidal position) of particle, modulo 2*pi.
c Phi is in R(2).

       phi = mod(r(2), tupi)
       if (phi .lt. 0.0d0) phi = phi + tupi
d       write (8, *) ' phi=', phi


c Now find the toroidal region in which particle lies
       itorreg = srchz (phi, phistt, ntorreg+1)
d       write (8, *) ' INTLIM: itorreg=', itorreg


       if ((itorreg .le. 0) .or. (itorreg .gt. ntorreg)) then
          print *, ' Particle toroidal angle (mod 2pi) is not in range:',
     .      phi, ' r(2)=', r(2), ' itorreg=', itorreg
          hitlim = .true.
          return
       endif

c Now find whether particle has hit limiter
c First, find interval in z along left limiter containing particle's z
       iz = srchz (r(3), zlimplf(1, itorreg), nlimplf(itorreg)) 
d      write (8,*) ' INTLIM: first search, iz, z, r=', iz, r(3), r(1)

       if (iz .le. 0) then
	   
c come here if particle's z coord is not in range of limiter values--
c means particle somehow left plasma region of vessel, so count it
c as hitting limiter
          print *, 'r(3) = ', r(3),
     >         ' outside inner z-range, tor. region = ', itorreg
         hitlim = .true.
         return
       endif

c If we've gotten here, it is because z is in range, and we can now 
c test the R of the particle by interpolating the line segment that
c describes the limiter outline here. Particle is
c between iz & iz+1.

       rint = (rlimplf(iz+1, itorreg) - rlimplf(iz, itorreg)) * 
     .   (r(3) - zlimplf(iz, itorreg)) / (zlimplf(iz+1, itorreg) 
     .   - zlimplf(iz, itorreg)) + rlimplf(iz, itorreg)

d      write(8,*) 'INTLIM: 1st RINT eval: ', rint,r(1)

       if (r(1) .le. rint) then
          print *, 'r(1) = ', r(1), 'interpolated : ', rint,
     >         ' outside inner r-range, tor. region = ', itorreg
          hitlim = .true.
          return
       endif

c Now repeat process for large-R limiter.

       iz = srchz (r(3), zlimprt(1, itorreg), nlimprt(itorreg))
d      write (8, *) ' INTLIM: 2nd search, iz= ', iz

       if (iz .le. 0) then
	   
c come here if particle's z coord is not in range of limiter values--
c means particle somehow left plasma region of vessel, so count
c it as hitting limiter
          print *, 'r(3) = ', r(3),
     >         ' outside outer z-range, tor. region = ', itorreg
         hitlim = .true.
         return
       endif

c If we've gotten here, it is because z is in range, and we can now 
c test the R of the particle by interpolating the line segment that
c describes the limiter outline here. Particle is
c between iz & iz+1.

       rint = (rlimprt(iz+1, itorreg) - rlimprt(iz, itorreg)) * 
     .   (r(3) - zlimprt(iz, itorreg)) / (zlimprt(iz+1, itorreg) 
     .   - zlimprt(iz, itorreg)) + rlimprt(iz, itorreg)

d      write (8,*) ' INTLIM: 2nd RINT eval: ', rint, r(1)

       if (r(1) .ge. rint) then
          print *, 'r(1) = ', r(1), 'interpolated : ', rint,
     >         ' outside outer r-range, tor. region = ', itorreg          
          hitlim = .true.
          return
       endif

       return
       end
