
      subroutine pmedinit
	  
c**********************************************************************
c**                                                                  **
c**     SUBPROGRAM DESCRIPTION:                                      **
c**             pmedinit calls weqdsk to initialize grid.
c**             now renamed rdeqdsk                                  **
c**                                                                  **
c**     RECORD OF MODIFICATION:                                      **
c**          03/07/89..........first created                         **
c**          03/07/89..........revised                               **
c**                                                                  **
c**                                                                  **
c**********************************************************************

      include 'dcalc.h'

      include 'eparmdu129.h'
      common/gtable/rgrid(nw),zgrid(nh)
      common/cpsi/psi(nwnh),psibry,simag,sidif,xpsi(nwnh),eouter
      common/efitcom/ifname, ifdir
      common/efitcom2/ifname2

      common/mw1com/mw1,mh1

      character*80 ifname, ifdir
      character*160 ifname2,filenm
      dimension workaa(nw)

      imfit1=1
      imfit2=1
c--------------------------------------------------------------------
c--  psi and B field calculations only                             --
c--------------------------------------------------------------------
 
c this makes ifname2 = gwell.dat where ifname is stored from the input file
c scanned by subroutine rdpar
	       
      print *, 'ifname: ', ifname
      print *, 'ifdir: ', ifdir
      ifname2=TRIM(ADJUSTL(ifdir))//'/g'//ifname
      print *, 'pmedinit: EQDSK file : '//ifname2
c      ifname2=ifname

      call rdeqdsk(ifname2,imfit1,mw1,mh1,ier)

      if (ier .ne. 0) print *, ' Return code from WEQDSK is nonzero:', ier

      return
      end
