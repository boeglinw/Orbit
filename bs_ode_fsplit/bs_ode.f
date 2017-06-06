c========================================================/bs_ode.for
c Bulirsch-Stoer Ordinary Differental Equation driver
c 26-Jan-94  Adapted BS_ODE code from NBretz Toray code.
c	Replaces IMSL V9 DVERK (Runge Kutta differential equation solver, 5th
c	and 6th order method).  See Numrical Recipes, ch 15.2

c
c Y is changed upon return
c S is changed upon return

      subroutine bs_ode(derv,		!user supplied derivative subroutine
     *			neqn,		!# of equations to be integrated
     *			y,		!array of dependent values
     *			s,		!independent variable,changed on return
     *			ds,		!stepsize
     *			relerr,		!error tolerances
     *		        iflag)		!0=ok,6=tolerance <=0.0


      implicit none

      save
      integer*4 neqn, iflag

      real*8 s, ds, relerr, omega
      real*8 y(neqn)
      external derv, rzextr, pzextr
      common/fcncom/ omega

c local variables
      integer*4 i, ncount, maxstep, n, idum, maxrk, j
      real*8 x, eps, htry, yscal(6), hdid, hnext, dydx(6), yy(6)
      real*8 totdid, x2, hmin, rdum
      data maxstep/ 25 /, maxrk/ 3 /

c -------------------- executable --------------------
c store current values
      x=s	!starting point

      htry=ds	!stepsize wanted to take

      x2=s+ds	!end step

      hmin=ds/maxstep	!minimum step size

      do i=1,neqn	!save original values
	yy(i)=y(i)
      end do
 


c try rational function extrapolation first, we want to go from s to s+ds

      ncount=0	!# of times hdid < htry
      totdid=0.	!total stepsizes taken


c set initial error
      iflag=0

c get eps
      eps=relerr
      if(eps.le.0.) then
	iflag=6
	write (6, *)' <BS_ODE> Error tolerance <= 0.0, unable to continue'
      end if

c 12/11/92  change to match Numerical Recipes ODEINT,pp.559-560

      do n=1,maxstep

c get derivatives, dydx
	call derv(idum,rdum,yy,dydx)

c set up error scaling vector

	do i=1,neqn
	  yscal(i)=abs(yy(i)) + abs(htry*dydx(i)) + 1.e-30
	end do

c check if we overshot end pt
	if( (x+htry-x2)*(x+htry-s) .gt.0.0) htry=x2-x

c try to take a step
    1	call bsstep(yy,dydx,neqn,x,htry,eps,yscal,hdid,hnext,derv,
     *	  rzextr)


	if( (x-x2)*(x2-s).ge.0.) then	!are we done?
	  goto 900	!-------------------->quit loop
	end if


c x is where we have moved to in hdid stepsize
c
c if bsstep failed to reach htry, then try a Runge-Kutta step
c starting at where bsstep left off

	if(hdid.lt.htry) then
ccjf	  if(i.eq.1) then	!first bsstep iteration?
	  ncount=ncount + 1
	  if(ncount.lt.10) then
	    htry=hnext	!use the next recommended B-S stepsize
	    goto 800	!loop again
	  else
	    htry=hnext/16.	!reduce step size for R-K
	    do j=1,maxrk
	      call derv(idum,rdum,yy,dydx)    !yy is changed, so need new dydx

	      do i=1,neqn
		yscal(i)=abs(yy(i)) + abs(htry*dydx(i))+1.e-30	!so no 0's
	      end do


c take a Runge_Kutta step
	      call rkqc(yy,dydx,neqn,x,htry,eps,yscal,hdid,hnext,derv)
	      htry=hnext
	    end do	!end runge kutta loop
	    ncount=0
	  end if	!if first bsstep loop or try runge kutta
	end if	!if didn't reach htry

c x and yy should be new values, hdid is how far we went
c now reset htry, loop back and try bsstep again

	htry=hnext
  800	continue
      end do	!end of maxsteps loop

c if we get here, we didn't converge

      write (6, *)' bs_ode: could not converge on step size ',ds,
     *  ' at this X location',x

c return new s, y
  900 s=x
      do i=1,neqn
	y(i)=yy(i)
      end do


      return
      end
