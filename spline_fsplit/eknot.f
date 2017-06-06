


      subroutine eknot(n,x,k,xk)
	  
c**********************************************************************
c**                                                                  **
c**                                                                  **
c**                                                                  **
c********************************************************************** 	  
c given the ordered data points x(1)<...<x(n), this subroutine generates
c a knot sequence with not-a-knot end conditions (like BSNAK from IMSL)
c Some of this is discussed in de Boor(1978), page 211.
C^d^  implicit integer*4 (i-n), real*8 (a-h, o-z)
      implicit real*8 (a-h, o-z)

      save
        dimension x(n),xk(n+k)
        integer kh
c
        do i=1,k
        xk(i)=x(1)
        ii=i+n
        xk(ii)= x(n)+1.e-5
        enddo
        kh=k/2
        k2=kh+kh
        if (k2.eq.k) then
c even k, place knots at data points
        do i=k+1,n
        xk(i)=x(i-kh)
        enddo
        else
c odd k, place knots in between data points
        do i=k+1,n
        xk(i)=.5*(x(i-kh)+x(i-1-kh))
        enddo
        end if
        return
        end
