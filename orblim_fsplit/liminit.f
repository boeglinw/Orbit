
       subroutine liminit
	   
********************************************************************
c This routine reads in the data which describes the non-axisymmetric
c limiter outline (as regions demarcated by toroidal angle).  This
c data is then used by the subroutine INTLIM which tests to see whether
c the particle has hit the limiter.
c*********************************************************************
   
      implicit real*8 (a-h, o-z)

      save
      common/limnam/limfnam, limfdir
      character*60 limfnam, limfdir
      parameter (nlimf = 17)

c     LIMFNAM is the name of the file which contains the limiter data
c     NLIMF is the unit number for reading the limiter data file

      parameter (ntorrmx = 30)
      parameter (nlplmx = 80)
      parameter (nlprmx = 80)

c     NTORRMX is the maximum number of toroidal regions (ie intervals in
c       phi, the toroidal angle).  Used to dimension arrays.
c     NLPLMX is the maximum number of points for defining the left part
c       (ie inboard part) of the limiter
c     NLPRMX maximum number of points for right (outboard) part of limiter

      common /lim/ ntorreg, phistt(ntorrmx+1), nlimplf(ntorrmx),
     .  nlimprt(ntorrmx), rlimplf(nlplmx, ntorrmx), 
     .  zlimplf(nlplmx, ntorrmx), rlimprt(nlprmx, ntorrmx),
     .  zlimprt(nlprmx, ntorrmx)
      integer ntorreg, nlimplf, nlimprt
      real*8 phistt, rlimplf, zlimplf, rlimprt, zlimptrt

c     NTORREG number of toroidal regions (ie intervals in phi, the toroidal
c       angle).  Must be .LE. NTORRMX
c     PHISTT toroidal angle at which each toroidal region starts (ending
c       angle is the same as the starting angle of the next region, and
c       there is one extra element at the end of this array, which is 
c       set, by definition, to 2*pi)
c     NLIMPLF, NLIMPRT number of points describing, respectively, the left
c       and right halves of limiter for this particular toroidal region.
c     RLIMPLF, ZLIMPLF R & z coordinates of points which define the limiter
c       on the left side (small R side) of the plasma for each toroidal 
c       region
c     RLIMPRT, ZLIMPRT same as above, but for right (large R) side of plasma

      character*4 cdum

c     CDUM for skipping comment lines in input file

      common /consts/ pm, echrg, emass, pi

      common /midlim/ rmplimo(2*ntorrmx+1), phmplimo(2*ntorrmx+1),
     .    rmplimi(2*ntorrmx+1), phmplimi(2*ntorrmx+1)
      real*4 rmplimo, phmplimo, rmplimi, phmplimi

c     RMPLIMO is the array of R values for the corners of the outer limiter at
c       the midplane
c     PHMPLIMO is the array of phi values (toroidal angles) for the corners
c       of the outer limiter at the midplane
c     RMPLIMI is the array of R values for the corners of the inner limiter at
c       the midplane
c     PHMPLIMI is the array of phi values (toroidal angles) for the corners
c       of the inner limiter at the midplane

      parameter (nidmax=10, nidmxm1=nidmax-1)
      common /limdrw/ nid, idrwseg(nidmax), ilintyp(nidmax)
      integer nid, idrwseg, ilintyp
      data nid /1/, idrwseg /1, nidmxm1*0/, ilintyp /0, nidmxm1*-1/

c     NIDMAX is the maximum number of limiter outlines to draw.
c     NID is the actual number of limiter outlines to draw.
c     IDRWSEG is the second subscript into RLIMPLF, ZLIMPLF, etc. of a
c       segment of the limiter which is to be drawn.  
c     ILINTYP is the SGLIB line type to use in the plot of this segment
c       (dashed, solid, dotted, etc.)

      parameter (tupi = 6.283185307180)

      integer srchz

      open (unit=nlimf, file=TRIM(ADJUSTL(limfdir))//'/'//limfnam, status='old')
      READ(nlimf, *) ntorreg

      print *, 'lim :', ntorreg

      if ((ntorreg .gt. ntorrmx) .or. (ntorreg .le. 0)) then
         print *, ' Number of toroidal regions for limiter, ', ntorreg,
     .     ' must be between 0 and ', ntorrmx, ', but is not.'
         stop
      endif
     
c     For number of toroidal regions specified, read starting toroidal angle, 
c     number of points (left hand side & right hand side) of data defining
c     limiter in this region, then read the actual coordinates.

      do i= 1, ntorreg
         READ(nlimf, '(a4)') cdum
c         print *, 'lim :', cdum
         READ(nlimf, *) phistt(i)
c         print *, 'lim :',i,phistt(i)
         if ((phistt(i) .lt. 0.0) .or. (phistt(i) .ge. 360.0)) then
            print *, ' Angle defining range of limiter is less than 0 or',
     .        ' greater than 360 degrees: ', phistt(i), ' (region ', i,
     .        ')'
            stop
         endif

c       Take input angles as degrees, but convert them now to radians.

         phistt(i) = phistt(i) * pi/180.
         if (i .eq. 1) then
            if (phistt(1) .ne. 0.0) then
               print *, ' Starting toroidal angle of first toroidal region',
     .              ' must be 0.0, but is not: ', phistt(1)*180./pi
               stop
            endif
         else
            if (phistt(i) .le. phistt(i-1)) then
               print *, ' Starting angles of toroidal regions must be ',
     .              ' in ascending order, but are not for region ', i
               stop
            endif
         endif

         READ(nlimf, *) nlimplf(i), nlimprt(i)

         print *, 'lim :',nlimplf(i), nlimprt(i)

         if (nlimplf(i) .gt. nlplmx) then
            print *, ' Number of limiter points on left, ', nlimplf(i),
     .        ' exceeds limit of ', nlplmx, ' for region ', i
            stop
         endif
         if (nlimprt(i) .gt. nlprmx) then
            print *, ' Number of limiter points on right, ', nlimprt(i),
     .        ' exceeds limit of ', nlprmx, ' for region ', i
            stop
         endif

         READ(nlimf, *) (rlimplf(j, i), j= 1, nlimplf(i))
         READ(nlimf, *) (zlimplf(j, i), j= 1, nlimplf(i))
         
         do j = 1,  nlimplf(i)
            print *, 'left-lim :', i, j, rlimplf(j, i), zlimplf(j, i)
         enddo

         READ(nlimf, *) (rlimprt(j, i), j= 1, nlimprt(i))
         READ(nlimf, *) (zlimprt(j, i), j= 1, nlimprt(i))

         do j = 1,  nlimprt(i)
            print *, 'rt-lim :', i, j, rlimprt(j, i), zlimprt(j, i)
         enddo

      enddo

      phistt(ntorreg+1) = tupi

c Now read info on how to display limiters in plot

      READ(nlimf, '(a4)') cdum
      READ(nlimf, *) nid
d     write (8, *) ' nid=', nid


      do i= 1, nid
         READ(nlimf, *) idrwseg(i), ilintyp(i)
d        write (8, *) ' idrwseg, ilintyp=', idrwseg(i), ilintyp(i)


      enddo

c That's all the limiter data, so close the file

      close (unit=nlimf)

d      write (8, *) ' ntorreg=', ntorreg
d      write (8, *) ' phistt=', (phistt(i), i= 1, ntorreg+1)
d      write (8, *) ' nlimplf=', (nlimplf(i), i= 1, ntorreg)
d      write (8, *) ' nlimprt=', (nlimprt(i), i= 1, ntorreg)
d      do i= 1, ntorreg
d         write (8, *) ' rlimplf=', (rlimplf(j, i), j= 1, nlimplf(i))
d         write (8, *) ' zlimplf=', (zlimplf(j, i), j= 1, nlimplf(i))
d         write (8, *) ' rlimplf=', (rlimprt(j, i), j= 1, nlimprt(i))
d         write (8, *) ' zlimplf=', (zlimprt(j, i), j= 1, nlimprt(i))
d      enddo

c Now find the R and phi values of the corners of the limiter at 
c the midplane

      do 250 itorreg= 1, ntorreg
d     write (8, *) ' In LIMINIT, midplane lim loop, itorreg=',itorreg

c       First do inner wall

        iz = srchz (0.0d0, zlimplf(1, itorreg), nlimplf(itorreg))

        rint = (rlimplf(iz+1, itorreg) - rlimplf(iz, itorreg)) * 
     .    (0.0d0 - zlimplf(iz, itorreg)) / (zlimplf(iz+1, itorreg) 
     .    - zlimplf(iz, itorreg)) + rlimplf(iz, itorreg)

d       write (8, *) ' iz, rint (left) =', iz, rint
        print *, ' iz, rint (left) =', iz, rint, 
     >       rlimplf(iz+1, itorreg), rlimplf(iz, itorreg)
        print *, ' iz,  (left) =', iz, 
     >       zlimplf(iz+1, itorreg), zlimplf(iz, itorreg)

        rmplimi(2*itorreg-1) = rint
        rmplimi(2*itorreg) = rint
        phmplimi(2*itorreg-1) = phistt(itorreg)
        phmplimi(2*itorreg) = phistt(itorreg+1)

c       Now do outer wall

        iz = srchz (0.0d0, zlimprt(1, itorreg), nlimprt(itorreg))
        rint = (rlimprt(iz+1, itorreg) - rlimprt(iz, itorreg)) * 
     .    (0.0d0 - zlimprt(iz, itorreg)) / (zlimprt(iz+1, itorreg) 
     .    - zlimprt(iz, itorreg)) + rlimprt(iz, itorreg)
d       write (8, *) ' iz, rint (right) =', iz, rint

        rmplimo(2*itorreg-1) = rint
        rmplimo(2*itorreg) = rint
        phmplimo(2*itorreg-1) = phistt(itorreg)
        phmplimo(2*itorreg) = phistt(itorreg+1)
		
250   continue

c     Now connect back at the end

      rmplimi(2*ntorreg+1) = rmplimi(1)
      rmplimo(2*ntorreg+1) = rmplimo(1)
      phmplimi(2*ntorreg+1) = 6.283185
      phmplimo(2*ntorreg+1) = 6.283185
      

d     write (8, *) ' In LIMINIT: ntorreg = ', ntorreg
d     do 290 i= 1, 2*ntorreg+1
d       write (8, *) ' (R, phi) for inner & outer limiter: ',
d    .    rmplimi(i), phmplimi(i), rmplimo(i), phmplimo(i)
d290  continue

c write the data for plotting
      open(52, file = 'limiter_drawing.data')
c R Z poins for the limiter regions 
      nmid = 2*ntorreg+1
      do i= 1, ntorreg
         write (52, *) '-region, ', i, phistt(i)
         write (52,*) '--left(r,z), ', nlimplf(i)
         do k = 1, nlimplf(i)
            write(52, *) rlimplf(k, i), zlimplf(k, i)
         enddo
         write (52,*) '--right(r,z), ', nlimprt(i)
         do k = 1, nlimprt(i)
            write(52, *) rlimprt(k, i), zlimprt(k, i)
         enddo
c midplane drawing
      enddo
      write (52,*) '--midplane, ', nmid
      write(52, *) '---inner(r,phi)'
      do k = 1, nmid
         write(52, *) rmplimi(k), phmplimi(k)
      enddo
      write(52, *) '---outer(r,phi)'
      do k = 1, nmid
         write(52, *) rmplimo(k), phmplimo(k)
      enddo
      close(52)

      return


      end
