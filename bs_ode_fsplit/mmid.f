
c%
c======================================================/mmid.for
c Modified midpoint step
c
c  See Numerical Recipes, ch 15.
c Dependent variable vector Y of length NVAR and its derivative vector
c DYDX are input at XS.  Also input is HTOT, the total step to be made, and
c NSTEP, the number of substeps to be used.  The output is returned as YOUT,
c which need not be a distinct array from Y;  if it is distinct, however,
c then Y and DYDX are returned undamaged.  DERIVS is a user supplied routine
c to find the derivatives.


      subroutine mmid(y,dydx,nvar,xs,htot,nstep,yout,derivs)

      implicit none

      save
      integer*4 nvar, nstep, istop
      real*8 y(nvar), dydx(nvar), yout(nvar)
      real*8 xs, htot
      external derivs


c local variables
      integer*4 nmax, i, n, idum
      parameter (nmax=6)
      real*8 ym(nmax), yn(nmax), h, x, h2, swap, rdum

c ---------------------- executable ------------------------

      h=htot/nstep	!stepsize this trip

      do i=1,nvar
	ym(i)=y(i)
	yn(i)=y(i)+h*dydx(i)	!first step
      end do

      x=xs+h
      call derivs(idum,rdum,yn,yout)	!use yout for temp storage of derivatives

      h2=2.*h

      do n=2,nstep	!general step
	do i=1,nvar
	  swap=ym(i)+h2*yout(i)
	  ym(i)=yn(i)
	  yn(i)=swap
	end do
	x=x+h
	call derivs(idum,rdum,yn,yout)

      end do	!end general step

      do i=1,nvar	!last step
	yout(i)=0.5*(ym(i)+yn(i)+h*yout(i))
      end do

      return
      end
