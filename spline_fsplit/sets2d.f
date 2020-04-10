


      subroutine sets2d(s,cs,x,nx,bkx,lx,y,ny,bky,ly,wk,ier)
	  
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
c      s     - nx by ny array containing the function values at (x,y).
c              This is a 1-d array, k=k=(i-1)*ny+j.
c
c      x, y  - (x,y) location, arrays of length nx and ny.
c
c  Outputs:
c
c      cs    - array of spline coefficients of dimension (kubicx,
c              lubicx,kubicy,lubicy).
c
c      bkx, bky - interval coefficients of length lubicx+1 and lubicy+1.
c
c      lx, ly -   number of terms in bkx and bky.
c
c      ier   - rror parameter.
c
c  Work arrays:
c
c      wk    - of dimension at least nx by ny.
c------------------------------------------------------------------------------
      implicit real*8 (a-h, o-z)

      save
      include 'eparmdu129.h'
c
c      dimension s(1), x(nx), y(ny), wk(nx,ny),
c     .          xknot(kubicx + nw), yknot(kubicy + nh),
c     .          cs(kubicx, lubicx, kubicy, lubicy),
c     .          bkx(lubicx + 1), bky(lubicy + 1)

      dimension s(1), x(nx), y(ny), wk(nx,ny),
     .          cs(kubicx, nx-kubicx+1, kubicy, ny-kubicy+1),
     .          bkx(nx-kubicx + 2), bky(ny-kubicy + 2),
     .          xknot(nxknot), yknot(nyknot)
c     .          xknot(kubicx + nx), yknot(kubicy + ny),
c
c  Set up knots:
c
      call eknot (nx, x, kubicx, xknot)		
      call eknot (ny, y, kubicy, yknot)			
c
c  Save the original, use the work array
c
      do 10 i=1,nx
      do 10 j=1,ny
         k=(i-1)*ny+j
  10     wk(i,j) = s(k)
c
c  Calculate spline coefficients:
c
      call spl2bc (x, y, nx,ny,xknot, yknot, wk)	
c
c  Coefficients stored in bkx, bky, and c:
c
      call spl2pp (nx,ny,xknot, yknot, wk, bkx, lx, bky, ly, cs)
c
      return
      end
