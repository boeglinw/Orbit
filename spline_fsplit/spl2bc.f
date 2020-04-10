


      subroutine spl2bc(rgrid,zgrid,mw,mh,rknot,zknot,copy)
	  
c**********************************************************************
c**                                                                  **
c**                                                                  **
c**                                                                  **
c********************************************************************** 	  
c calculates the b-spline coeficients
C^d^  implicit integer*4 (i-n), real*8 (a-h, o-z)

      implicit real*8 (a-h, o-z)

      save
      parameter (nw=129,nh=129,krord=4,kzord=4)
      dimension rgrid(mw),zgrid(mh)
c      dimension rknot(nw+krord),zknot(nh+kzord),copy(mw,mh)
      dimension rknot(mw+krord),zknot(mh+kzord),copy(mw,mh)

c------------------------------------------------------------------
c-- change dimension of work2 and work3 from nw to nh            --
c-- to ensure the cases when nh > nw     ll, 93/04/01            --
c------------------------------------------------------------------
      
c      dimension work1(mw,mh),work2(mh),work3(mh,2*krord-1)
      dimension work1(nw,nh),work2(nh),work3(nh,2*krord-1)

      call spli2d(rgrid,copy,rknot,mw,krord,mh,work2,work3,work1,iflag)
      if (iflag.ne.1) print*,' error in first spli2d, iflag=',iflag

      call spli2d(zgrid,work1,zknot,mh,kzord,mw,work2,work3,copy,iflag)
      if (iflag.ne.1) print*,' error in second spli2d, iflag=',iflag
      return
      end
