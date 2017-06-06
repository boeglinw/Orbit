      subroutine seva2d(bkx,lx,bky,ly,cs,nx,ny,xl,yl,fs,ier,icalc)

c**********************************************************************
c**                                                                  **
c**                                                                  **
c**                                                                  **
c********************************************************************** 	  
c------------------------------------------------------------------------------
c--  S.Thompson  92/05/18                                                    --
c--    Bicubic spline routines.                                              --
c--    Put together with routines from E.Solano.                             --
c--  SMWolfe     93/12/17                                                    --
c--    Modifed to avoid over-writing the original array.                     --
c--              94/04/28                                                    --
c--    Updated.                                                              --
c------------------------------------------------------------------------------
c  Inputs:
c
c      cs       - array of spline coefficients of dimension (kubicx,
c                 lubicx,kubicy,lubicy) from sets2d.
c
c      bkx, bky - interval coefficients of length lubicx+1 and lubicy+1 from
c                 sets2d.
c
c      lx, ly   - number of terms in bkx and bky from sets2d.
c
c      xl, yl   - the point at which interpolations are desired.
c
c      nx, ny   - grid dimensions
c
c  Outputs:
c
c      fs       - vector containing results depending on icalc:
c                 icalc              fs
c                   1                f
c                   2                fx
c                   3                fy
c                   4                fxy
c                   5                fxx
c                   6                fyy
c
c      ier      - error parameter.
c
c-------------------------------------------------------------------------------
      implicit real*8 (a-h, o-z)
      
      save
      include 'eparmdu129.h'
c
      integer ier, lx, ly
c      real cs(kubicx,lubicx,kubicy,lubicy),xl,yl,fs(6),bkx(1),bky(1)
      real*8 xl,yl,fs(6),bkx(1),bky(1)
      dimension cs(kubicx,nx-kubicx+1,kubicy,ny-kubicy+1) 
c
c  Local Variable Specifications:
c
      dimension work0(4),work1(4),work2(4)
      data n00/0/,n11/1/,n22/2/
c	  integer n00, n11, n22
c
c  Evaluate function and its partial derivatives at (XL, YL):
c
c  First do all the lookup and interpolation stuff.
c  This is the most time consuming part of the evaluation, so
c  don't do more than needed.
c
      call interv(bky,ly,yl,lef,mflag)
      call interv(bkx,lx,xl,ibk,ndummy)
      h = xl - bkx(ibk)
      do 41 jj=1,4
         work0(jj) = ppvalw(cs(1,ibk,jj,lef),h,n00)
         if (icalc.eq.1) goto 41
         work1(jj) = ppvalw(cs(1,ibk,jj,lef),h,n11)
         if (icalc.le.4) goto 41
         work2(jj) = ppvalw(cs(1,ibk,jj,lef),h,n22)
 41   continue
      h = yl - bky(lef)
      fs(1) = ppvalw(work0,h,n00)
      if (icalc.eq.1) return
      fs(2) = ppvalw(work1,h,n00)
      if (icalc.eq.2) return
      fs(3) = ppvalw(work0,h,n11)
      if (icalc.eq.3) return
      fs(4) = ppvalw(work1,h,n11)
      if (icalc.eq.4) return
      fs(5) = ppvalw(work2,h,n00)
      if (icalc.eq.5) return
      fs(6) = ppvalw(work0,h,n22)
C
      return
      end
