	  
	  
	  
      function ppvalw (coef, x, jd )

c**********************************************************************
c**                                                                  **
c**                                                                  **
c**                                                                  **
c********************************************************************** 
C-----------------------------------------------------------------------        
C  Modified for optimization by S.J. Thompson, 30-Aug-1993
c  Revised to eliminate call to interv by S.M.Wolfe, 17-Dec-1993
c          and to use ASF's for evaluation
c  This routine performs only the innermost guts of the spline evaluation
c  Assumes k=4 (cubic spline only). No other cases considered.
c does not call  interv
calculates value at  x  of  jd-th derivative of pp fct from pp-repr
c
c******  i n p u t  ****** to PPVALU, on which this is based.
c  break, coef, l, k.....forms the pp-representation of the function  f
c        to be evaluated. specifically, the j-th derivative of  f  is
c        given by
c
c     (d**j)f(x) = coef(j+1,i) + h*(coef(j+2,i) + h*( ... (coef(k-1,i) +
c                             + h*coef(k,i)/(k-j-1))/(k-j-2) ... )/2)/1
c
c        with  h = x - break(i),  and
c
c       i  =  max( 1 , max( j ;  break(j) .le. x , 1 .le. j .le. l ) ).
c
c  x.....the point at which to evaluate.
c        as used here, x is the distance from the break, not the absolute 
c        position. 
c  jd.....integer*4 giving the order of the derivative to be evaluat-
c        ed.  a s s u m e d  to be zero or positive.
c
c******  o u t p u t  ******
c  ppvalw.....the value of the (jd)-th derivative of  f  at  x.
c
c******  m e t h o d  ******
c     the interval index  i , appropriate for  x , is found through a
c  call to  interv . the formula above for the  jd-th derivative
c  of  f  is then evaluated (by nested multipication).
c
C-----------------------------------------------------------------------        
C   Variable declarations.
C-----------------------------------------------------------------------        
C^d^  implicit integer*4 (i-n), real*8 (a-h, o-z)
      implicit real*8 (a-h, o-z)

      save
      real*8 ppvalw,x
      dimension coef(4)
c----------------------------------------------------------------------
c ASF's may be slightly more efficient than the alternative
c----------------------------------------------------------------------
      d2(xx) = coef(4)*xx + coef(3)
      d1(xx) = (coef(4)*xx/2. + coef(3))*xx + coef(2)
      d0(xx) = ((coef(4)*xx/3. + coef(3))*xx/2. + 
     .           coef(2))*xx + coef(1)
C-----------------------------------------------------------------------        
C   Derivatives of order k or higher are identically zero.
C-----------------------------------------------------------------------        
C   Evaluate jd-th derivative of i-th polynomial piece at x .
C-----------------------------------------------------------------------        
      goto (1,2,3) jd+1
      ppvalw = 0.
      print *, 'Error (ppvalw): JD must be 0, 1, or 2.'
      print *, 'Execution terminated.'
      return
 1    ppvalw = d0(x)	! k = 4 , jd = 0
      return
 2    ppvalw = d1(x)	! k = 4 , jd = 1
      return
 3    ppvalw = d2(x)	! k = 4 , jd = 2
      return
      end
