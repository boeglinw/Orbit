	  
	  
	  
      subroutine banfac ( a, nrow, n, ndiag, middle, iflag )

c**********************************************************************
c**                                                                  **
c**                                                                  **
c**                                                                  **
c********************************************************************** 
C^d^  implicit integer*4 (i-n), real*8 (a-h, o-z)
      implicit real*8 (a-h, o-z)

      save
      dimension a(nrow,ndiag)
      iflag = 1
      ilo = middle - 1
      if (ilo)                          999,10,19
   10 do 11 i=1,n
         if(a(i,1) .eq. 0.)             go to 999
   11    continue
                                        return
   19 ihi = ndiag - middle
      if (ihi)                          999,20,29
   20 do 25 i=1,n
         if (a(i,middle) .eq. 0.)       go to 999
         jmax = min(ilo,n-i)
         if (jmax .lt. 1)               go to 25
         do 23 j=1,jmax
   23       a(i+j,middle-j) = a(i+j,middle-j)/a(i,middle)
   25    continue
                                        return
   29 do 50 i=1,n
         diag = a(i,middle)
         if (diag .eq. 0.)              go to 999
         jmax = min(ilo,n-i)
         if(jmax .lt. 1)                go to 50
         kmax = min(ihi,n-i)
         do 33 j=1,jmax
            mmj = middle-j
            a(i+j,mmj) = a(i+j,mmj)/diag
            do 33 k=1,kmax
   33          a(i+j,mmj+k) = a(i+j,mmj+k) - a(i+j,mmj)*a(i,middle+k)
   50    continue
                                        return
  999 iflag = 2
                                        return
      end
