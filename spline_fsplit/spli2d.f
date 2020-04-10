		
		
		
      subroutine spli2d ( tau, gtau, t, n, k, m, work, q, bcoef, iflag )
	  
c**********************************************************************
c**                                                                  **
c**                                                                  **
c**                                                                  **
c********************************************************************** 
  
c  calls bsplvb, banfac/slv
c  this is an extended version of  splint , for the use in tensor prod-
c  uct interpolation.
c
c   spli2d  produces the b-spline coeff.s  bcoef(j,.)  of the spline of
c   order  k  with knots  t (i), i=1,..., n + k , which takes on the
c   value  gtau (i,j)  at  tau (i), i=1,..., n ; j=1,..., m .
c
c******  i n p u t  ******
c  tau   array of length  n , containing data point abscissae.
c  a s s u m p t i o n . . .  tau  is strictly increasing
c  gtau(.,j)  corresponding array of length  n , containing data point
c        ordinates, j=1,...,m
c  t     knot sequence, of length  n+k
c  n     number of data points and dimension of spline space  s(k,t)
c  k     order of spline
c  m     number of data sets
c
c******  w o r k   a r e a  ******
c  work  a vector of length  n
c
c******  o u t p u t  ******
c  q     array of order  (n,2*k-1), containing the triangular factoriz-
c        ation of the coefficient matrix of the linear system for the b-
c        coefficients of the spline interpolant.
c           the b-coeffs for the interpolant of an additional data set
c        (tau(i),htau(i)), i=1,...,n  with the same data abscissae can
c        be obtained without going through all the calculations in this
c        routine, simply by loading  htau  into  bcoef  and then execut-
c        ing the    call banslv ( q, n, n, 2*k-1, k, bcoef )
c  bcoef the b-coefficients of the interpolant, of length  n
c  iflag an integer indicating success (= 1)  or failure (= 2)
c        the linear system to be solved is (theoretically) invertible if
c        and only if
c              t(i) .lt. tau(i) .lt. tau(i+k),    all i.
c        violation of this condition is certain to lead to  iflag = 2 .
c
c******  m e t h o d  ******
c     the i-th equation of the linear system  a*bcoef = b  for the b-co-
c  effs of the interpolant enforces interpolation at  tau(i), i=1,...,n.
c  hence,  b(i) = gtau(i), all i, and  a  is a band matrix with  2k-1
c  bands (if it is invertible).
c     the matrix  a  is generated row by row and stored, diagonal by di-
c  agonal, in the  c o l u m n s  of the array  q , with the main diag-
c  onal going into column  k .  see comments in the program below.
c     the banded system is then solved by a call to  banfac (which con-
c  structs the triangular factorization for  a  and stores it again in
c   q ), followed by a call to  banslv (which then obtains the solution
c   bcoef  by substitution).
c     banfac  does no pivoting, since the total positivity of the matrix
c  a  makes this unnecessary.
c
c      integer iflag,k,m,n,i,ilp1mx,j,jj,kpkm1,left,np1
c      real bcoef(m,n),gtau(n,m),q(n,7),t(n+k),tau(n),work(n),taui
C^d^  implicit integer*4 (i-n), real*8 (a-h, o-z)
      implicit real*8 (a-h, o-z)

      save
      dimension bcoef(m,n),gtau(n,m),q(n,2*k-1),t(n+k),tau(n),work(n)
c
      nnn=1
      np1 = n + 1
      kpkm1 = 2*k - 1
      left = k
c
c  ***   loop over  i  to construct the  n  interpolation equations
      do 30 i=1,n
         iindex=i
         taui = tau(iindex)
         ilp1mx = min(iindex+k,np1)
c        *** zero out all entries in row  i  of  a (in the 2k-1 bands)
         do 13 j=1,kpkm1
   13       q(iindex,j) = 0.
c        *** find  left  in the closed interval (i,i+k-1) such that
c                t(left) .le. tau(i) .lt. t(left+1)
c        matrix is singular if this is not possible
         left = max(left,i)
         if (taui .lt. t(left))         go to 998
   15       if (taui .lt. t(left+1))    go to 16
            left = left + 1
            if (left .lt. ilp1mx)       go to 15
         left = left - 1
         if (taui .gt. t(left+1))       go to 998
c        *** the i-th equation enforces interpolation at taui, hence
c        a(i,j) = b(j,k,t)(taui), all j. only the  k  entries with  j =
c        left-k+1,...,left actually might be nonzero. these  k  numbers
c        are returned, in  work  (used for temp.storage here), by the
c        following
   16    call bsplvb ( t, k, nnn, taui, left, work )
c        we therefore want  work(j) = b(left-k+j)(taui) to go into
c        a(i,left-k+j), i.e., into  q(i,left-i+j), since the i-th row of
c        a  is so stored in the i-th row of  q  that the (i,i)-entry of
c        a  goes into the  k-th  entry of  q.
         jj = left - iindex
         do 29 j=1,k
            jj = jj+1
            q(iindex,jj) = work(j)
   29    continue
   30    continue
c
c     ***obtain factorization of  a  , stored again in  q.
      call banfac ( q, n, n, kpkm1, k, iflag )
                                        go to (40,999), iflag
c     *** solve  a*bcoef = gtau  by backsubstitution
   40 do 50 j=1,m
         do 41 i=1,n
   41       work(i) = gtau(i,j)
         call banslv ( q, n, n, kpkm1, k, work )
         do 50 i=1,n
   50    bcoef(j,i) = work(i)
                                        return
  998 iflag = 2
  999 print 699
  699 format(41h linear system in  splint  not invertible)
                                        return
      end
