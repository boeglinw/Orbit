
c%
c=========================================================/rkqc
c Fifth order Runge-Kutta step with monitoring of local truncation error
c to ensure accuracy and adjust stepsize.  
c From Numerical Recipes, ch 15.

c 29-Mar-90 JFelt

c Input:
c	Y	- dependent variable vector
c	N	- size of vector
c	DYDX	- derivatives of Y
c	X	- independent variable
c	HTRY	- step size to be attempted (scaled down by 16 or 32 as nec.)
c	EPS	- required accuracy
c	YSCAL	- vector of error scaling factors (|y(i)|+|h*dydx(i)|+1.e-30)
c	DERIVS 	- user supplied subroutine to compute right-hand side derivs
c
c Output:
c	Y,X	- replaced with new values
c	HDID	- stepsize actually accomplished
c	HNEXT	- estimated next stepsize


      subroutine rkqc(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)

      implicit none

      save
      integer*4 n, istop
      real*8 x, htry, eps, hdid, hnext
      real*8 y(n), yscal(n), dydx(n)
      external derivs


c local variables
      integer*4 nmax, i, idum
      real*8 pgrow, pshrnk, fcor, one, safety, errcon, errmax, h, hh
      parameter ( nmax=6, pgrow=-0.20, pshrnk=-0.25, fcor=1./15.,
     *	one=1., safety=0.9, errcon=6.e-4)
      real*8 xsav, ytemp(nmax), ysav(nmax), dysav(nmax), rdum

c note:  errcon=(4/safety)**(1/pgrow)

c-----------------------executable-----------------------

c save initial values

      xsav=x
      do i=1,n
	ysav(i)=y(i)
	dysav(i)=dydx(i)
      end do

c set stepsize to initial trial value
      h=htry

c take two half steps

    1 hh=0.5*h
      call rk4(ysav,dysav,n,xsav,hh,ytemp,derivs)

      x=xsav+hh
      call derivs(idum,rdum,ytemp,dydx)


      call rk4(ytemp,dydx,n,x,hh,y,derivs)

      x=xsav+h

      if(x.eq.xsav) then
	write (6, *)' <rkqc> stepsize not significant'
      end if

c take the large step
      call rk4(ysav,dysav,n,xsav,h,ytemp,derivs)

c evaluate accuracy
      errmax=0.
      do i=1,n
	ytemp(i)=y(i)-ytemp(i)	! error estimate
	errmax=max(errmax,abs(ytemp(i)/yscal(i)))
      end do

      errmax=errmax/eps	!scale relative to required tolerance
      if(errmax.gt.one) then	!truncation error too large, reduce stepsize
	h=safety*h*(errmax**pshrnk)
	goto 1	!try again
      else	!step succeeded, compute next stepsize
	hdid=h
	if(errmax.gt.errcon) then
	  hnext=safety*h*(errmax**pgrow)
	else
	  hnext=4.*h
	end if
      end if

c mop up fifth-order truncation error
      do i=1,n
	y(i)=y(i)+ytemp(i)*fcor
      end do

      return
      end

