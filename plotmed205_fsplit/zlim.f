
      subroutine zlim(zero,nw,nh,limitr,xlim,ylim,x,y)
c**********************************************************************
c**                                                                  **
c**     MAIN PROGRAM:  MHD FITTING CODE                              **
c**                                                                  **
c**                                                                  **
c**     SUBPROGRAM DESCRIPTION:                                      **
c**          zlim determines whether points on the (x,y) grid are    **
c**          inside or outside of the boundary set by the limiters.  **
c**                                                                  **
c**                                                                  **
c**     CALLING ARGUMENTS:                                           **
c**       zero............1 if inside and 0 otherwise (returned)     **
c**       nw..............dimension of x (input)                     **
c**       nh..............dimension of y (input)                     **
c**       limitr..........number of limiter points (input)           **
c**       xlim............r coordinates of limiter (input)           **
c**       ylim............z coordinates of limiter (input)           **
c**       x...............r grid (input)                             **
c**       y...............z grid (input)                             **
c**                                                                  **
c**     REFERENCES:                                                  **
c**          (1)                                                     **
c**          (2)                                                     **
c**                                                                  **
c**     RECORD OF MODIFICATION:                                      **
c**          26/04/83..........first created                         **
c**                                                                  **
c**                                                                  **
c**                                                                  **
c**********************************************************************

      implicit real*8 (a-h, o-z)

      dimension  zero(1),xlim(1),ylim(1),x(1),y(1)
      logical first
      data first/.true./

      if (first) then
         open( 51, file = 'flux_limit.data')
         write (51, *) '# xlim : r - coordinate  of limiter'
         write (51, *) '# ylim : z - coordinate  of limiter'
         write(51, *) '#! xlim[f,0]/ ylim[f,1]/ '
         do i = 1, limitr
            write (51,*) xlim(i), ylim(i)
         enddo
         close(51)
         first = .false.
      endif
      kk = 0
      do 100 i = 1,nw
        do 100 j = 1,nh
          kk = kk + 1
          zero(kk) = 1.
          ncross = 0
          do 20 k = 1,limitr-1
            if ((ylim(k).lt.y(j)) .and. (ylim(k+1).lt.y(j))) go to 20
            if (x(i) .eq. xlim(k))  go to 20
              t = x(i) - xlim(k)
              s = xlim(k+1) - x(i)
              if ((t*s) .lt. 0.) go to 20
                di = (ylim(k+1)-ylim(k)) / (xlim(k+1)-xlim(k))
                f = ylim(k) + di*(x(i)-xlim(k))
                if (f .lt. y(j)) go to 20
                  ncross = ncross + 1
   20     continue
          mcross = .5*ncross
          mcross = 2*mcross
          if (ncross .eq. mcross) zero(kk) = 0.
  100 continue
      return
      end
