
c%
c======================================================/bsstep.for
c Roland Bulirsch and Josef Stoer Ordinary Differential Equations Solver
c See Numerical Recipes, ch. 15, Integration of Ordinary Differential Equations

c The following is an excerpt.
c Bulirsch-Stoer method:
c	1.  Consider the final answer of a numerical calculation as itself 
c	being an analytic function of an adjustable parameter like the stepsize
c	"h".  That analytic function can be probed by performing the 
c       calculation with various values of h, none of them being necessarily 
c	small enough to yield the accuracy desired.  The function is then 
c	fitted to some analytic form and evaluated at h=0.
c	2.  Rational function extrapolation in Richardson-type applications
c	used to  break the shackles of the power series and its limited radius
c	of convergence, out only to the distance of the first pole in the
c	complex plane.  Rational function fits can remain good approximations
c	to analytic functions even after the various terms in powers of h
c	all have comparable magnitudes.  In other words, h can be so large as
c	to make the whole notion of the "order" of the method meaningless--
c	and the method can still work superbly.
c	3. Use a method whose error function is strictly even, allowing the
c	rational function approximation to be in terms of the variable h**2
c	instead of just h.
c
c	A single Bulirsch-Stoer step takes us from x to x+H, where H is 
c	supposed to be quite a large (not at all infinitesimal) distance.  That
c	single step is a grand leap consisting of many substeps of modified
c	midpoint method, which are then rational-function extrapolated to zero
c	stepsize.

c	The sequence of separate attempts to cross the interval H is made with
c	increasing values of n, the number of substeps.  A conventional
c	sequence of n's is n=2,4,6,8,12,16,24,32,48,64,96,...,n(i)=2n(i-2)
c	(The more obvious choice n=2,4,8,16,... makes h too small too rapidly.)
c	For each step, we do not know in advance how far up this sequence we
c	will go.  After each successive n is tried, a rational function extra-
c	polation is attempted.  That extrapolation returns both extrapolated 
c	values and error estimates.  If the errors are not satisfactory, we go
c	higher in n.  If they are satisfactory, we go on to the next step and
c	begin anew with n=2.

c	Another adjustable parameter is the number of previous estimates of
c	the functions at x+H (different values of n) to incorporate into the
c	rational function fit.  One possibility would be to use all the
c	estimates available, n=2,4,6,...etc.  By the time n gets moderately
c	large, the early values are not very relevant, so use NUSE=7, the
c	the number of estimates to use as suggested by Gear.
c	Error control is enforced by monitoring internal consistency, and
c	adapting stepsize to match a prescribed bound on the local truncation
c	error.  Each new result from the sequence of modified midpoint integra-
c	tions allow an extension by one additional set of diagonals.  The size 
c	of the new correction added at each stage is taken as the conservative
c	error estimate.

c	Each Bulirsch-Stoer step effectively ranges of a factor of,say, up to
c	96/2 in stepsize.  Furthermore, each step can be effectively of very
c	high order in h**2--if "order" means anything at all for values of h
c	as large as the method often takes.  The range of accuracy already
c	accessible to a step, 48**2 to some high power, is immense.  It is a 
c	fairly rare event for a B-S step to fail at all, and such failure is 
c	usually associated with starting or ending transients or internal
c	singularities in the solution.

c	One desideratum is to keep the number of sequences tried below or at
c	the value NUSE, since after this point early information is thrown
c	away without affecting the solution.  As each step succeeds, we learn
c	the value of n at which success occurs.  We can then try the next 
c	step with H scaled so that the same value h will be reached on, say,
c	sequence number NUSE-1.  When NUSE-1 succeeds, we might tweak the step-
c	size up a little bit; while when NUSE succeeds, we might tweak it down
c	a bit.  If the step truly fails--goes all the way upto IMAX without
c	finding an acceptable error, then H must be decreased substantially and
c	the step retried.

c 	Keep a watch for failed steps, where HDID is returned with a value
c	less than HTRY.  If these are not rare, then Bulirsch-Stoer is in
c	trouble.  Can try to save it with either of two options:
c	- Take a couple of quality-controlled Runge_Kutta steps to get over the
c	  rough spot.  A call to RKQC is exactly substitutable for a call to
c	  BSSTEP.  (You might, however, reduce the suggested stepsize HTRY by
c	  a factor of 16 or 32, in recognition of Runge-Kutta's necessarily
c	  smaller steps.  If this isn't done, RKQC will have to seek out
c	  the smaller stepsize by itself, with additional effort.)
c	- If the problem is not a rough spot in the solutions, but purely an
c	  artifact in the rational function extrapolation, then a less drastic
c	  therapy is to use polynomial extrapolation instead of rational
c	  function extrapolation for a step or two.  Polynomial extrapolation
c	  does not involve any divisions, and it can be less finicky than
c	  rational extrapolation--also less powerful!


c
c Bulirsch-Stoer step with monitoring of local truncation error to ensure
c accuracy and adjust stepsize.  
c	Input :
c		Y, dependent variable vector of length, nv
c		DYDX, derivative at the starting value of the independent
c			variable, X
c		X, independent variable
c		HTRY, stepsize to be attempted
c		EPS, required accuracy
c		YSCAL, vector against which the error is scaled
c		DERIVS, user supplied subroutine to compute derivatives
c		EXTRAP, user supplied extrapolation method
c
c	Output:
c		Y,X are replaced with new values
c		HDID, is the stepsize which was actually accomplished
c		HNEXT, estimated next stepsize


c 20-Jun-90  add common block stp to pass back istop flag from derivs.

      SUBROUTINE BSSTEP(Y,DYDX,NV,X,HTRY,EPS,YSCAL,HDID,HNEXT,DERIVS,
     *		EXTRAP)

      implicit none
      save
      integer*4 nv
      real*8 y(nv), x, htry, eps, yscal(nv), hdid, hnext, dydx(nv)
      external derivs, extrap


c local variables
      integer*4  nmax, imax, nuse
      real*8 one, shrink, grow, xsav, xest, errmax, h
      PARAMETER (NMAX=10,IMAX=11,NUSE=7,ONE=1.E0,SHRINK=.95E0,GROW=1.2E0
     *)

      integer*4 nseq(imax), i, j
      real*8 yerr(nmax), ysav(nmax), dysav(nmax), yseq(nmax)


      DATA NSEQ /2,4,6,8,12,16,24,32,48,64,96/

c ----------------------- executable -----------------

c save starting values

      H=HTRY
      XSAV=X
      DO 11 I=1,NV
        YSAV(I)=Y(I)
        DYSAV(I)=DYDX(I)
11    CONTINUE

c evaluate the sequence of modified midpoint integrations

1     DO 10 I=1,IMAX
        CALL MMID(YSAV,DYSAV,NV,XSAV,H,NSEQ(I),YSEQ,DERIVS)

        XEST=(H/NSEQ(I))**2	!squared since error series is even

c perform rational function extrapolation or polynomial extrapolation
c new y's returned

        CALL extrap(I,XEST,YSEQ,Y,YERR,NV,NUSE)


        ERRMAX=0.	!check local truncation error
        DO 12 J=1,NV
          ERRMAX=MAX(ERRMAX,ABS(YERR(J)/YSCAL(J)))
12      CONTINUE

c scale accuracy relative to tolerance
        ERRMAX=ERRMAX/EPS
        IF(ERRMAX.LT.ONE) THEN		!step converged
          X=X+H		!new x
          HDID=H
          IF(I.EQ.NUSE)THEN
            HNEXT=H*SHRINK
          ELSE IF(I.EQ.NUSE-1)THEN
            HNEXT=H*GROW
          ELSE
            HNEXT=(H*NSEQ(NUSE-1))/NSEQ(I)
          ENDIF
          RETURN	!normal return
        ENDIF
10    end do	!end imax loop

c if here, then step failed, quite unusual for this method
c reduce stepsize and try again

      H=0.25*H/2**((IMAX-NUSE)/2)

      if(x+h .eq. x) then
	write (6, *)
     *	' <BSSTEP> Step size underflow, retrying with smaller stepsize'

      end if
      GOTO 1
      END
