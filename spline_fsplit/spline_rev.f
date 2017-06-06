
c
c   This routine is required if the CVS revision numbers are to 
c   survive an optimization.
c
c
c   $Date: 1997/04/05 01:43:16 $ $Author: peng $
c
 
	     
	  subroutine spline_rev(i)
	  
c**********************************************************************
c**                                                                  **
c**                                                                  **
c**                                                                  **
c********************************************************************** 	  
      CHARACTER*100 opt
      character*10 s 
      if( i .eq. 0) s = 
     .'@(#)$RCSfile: spline.f,v $ $Revision: 2.1 $\000'
      return
      end
