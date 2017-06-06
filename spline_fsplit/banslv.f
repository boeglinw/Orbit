


      subroutine banslv ( a, nrow, n, ndiag, middle, b )
	  
c**********************************************************************
c**                                                                  **
c**                                                                  **
c**                                                                  **
c********************************************************************** 	  
C^d^  implicit integer*4 (i-n), real*8 (a-h, o-z)
      implicit real*8 (a-h, o-z)

      save
      dimension a(nrow,ndiag),b(n)
      if (n .eq. 1)                     go to 21
      ilo = middle - 1
      if (ilo .lt. 1)                   go to 21
      do 19 i=2,n
         jmax = min(i-1,ilo)
         do 19 j=1,jmax
   19       b(i) = b(i) - b(i-j)*a(i,middle-j)
c
   21 ihi = ndiag-middle
      do 30 i=n,1,-1
         jmax = min(n-i,ihi)
         if (jmax .lt. 1)               go to 30
         do 25 j=1,jmax
   25       b(i) = b(i) - b(i+j)*a(i,middle+j)
   30    b(i) = b(i)/a(i,middle)
                                        return
      end
