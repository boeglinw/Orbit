
c%
c=======================================================/rk4
c Given values for N variables Y and their derivatives DYDX known at X,
c use the fourth-order Runge-Kutta method to advance the solution over an
c interval H and return the incremented variables as YOUT, which need not
c be a distinct array from Y.  The user supplies the subroutine 
c DERIVS(n,x,y,dydx) which returns derivatives DYDX at X.

c From Numerical Recipes, ch 15

      subroutine rk4(y,dydx,n,x,h,yout,derivs)
      implicit none

      save
      integer*4 n
      real*8 x, h, y(n), dydx(n), yout(n)
      external derivs


c local variables
      integer*4 nmax, i, idum
      parameter (nmax=6)
      real*8 hh, h6, xh, rdum
      real*8 yt(nmax), dyt(nmax), dym(nmax)

c---------------------- executable ------------------

      hh=h*0.5
      h6=h/6.
      xh=x+hh

c first step
      do i=1,n
	yt(i)=y(i)+hh*dydx(i)
      end do

c second step
      call derivs(idum,rdum,yt,dyt)

      do i=1,n
	yt(i)=y(i)+hh*dyt(i)
      end do

c third step
      call derivs(idum,rdum,yt,dym)

      do i=1,n
	yt(i)=y(i)+h*dym(i)
	dym(i)=dyt(i)+dym(i)
      end do

c fourth step
      call derivs(idum,rdum,yt,dyt)


c accumulate increments with proper weights

      do i=1,n
	yout(i)=y(i) + h6*(dydx(i) + dyt(i) + 2.*dym(i))
      end do

      return
      end
