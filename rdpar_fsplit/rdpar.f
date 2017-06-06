       SUBROUTINE rdpar
C*********************************************************************
c
c    This subroutine reads all the parameters from the input file.
c  (This stuff used to be in the main program.)
c
C*********************************************************************

      include 'dcalc.h'
      save

      common/efitcom/ifname
      common/limnam/limfnam
      character*60 limfnam
      common/magcom/icur,cur,bt
      common/geocom/rmax,V0,W1,W2
      common/p14com/zdet,adet,ener 
 
      parameter (mxorpt=64000)
      common/data5/x(mxorpt+1),y(mxorpt+1),z(mxorpt+1),
     .    x2(mxorpt+1),y2(mxorpt+1)
      common/switch/ifor,ipl1,ipl2,irmn,rarray(201),rmn
      common/char/fname,coment
      common /versn/ codenm, compdt, comptm
      common/iter/iters(5),start(5),stop(5),step(5),name(5)
      dimension save(5)
      common /data6/ iflxplt, steppsi(mxorpt+1), pphi(mxorpt+1),
     .    mu(mxorpt+1), energy(mxorpt+1), b_r(mxorpt+1), b_phi(mxorpt+1), 
     .    b_z(mxorpt+1)
      real*8 mu
      common/shftcm/is

      logical eq99
      real*8 is
      character*6 ifil
      character*80 fname
      character*80 ifname
      character*80 ifdir
      character*80 limfdir 
      character*80 coment
      character*9  compdt, comptm
      character*35 codenm
      character*9 rundt, runtm
c
c     coment      comment field describing this run of code
c     codenm      name of this version of code
c     compdt      compilation date of code (manual for now, auto later?)
c     comptm      compilation time of code (")
c     rundt       date of current run of code
c     runtm       time of current run of code
c

      common /magsrc/ magfil
      character*60 magfil

c     MAGFIL	  name of file containing CHS coil currents

c                        The following common block is for 2 plots per page!

      common /plot2/ ixlpol2, ixhpol2, jybpol2, jytpol2, ixlpln2, ixhpln2,
     .   jybpln2, jytpln2

c  IXLPOL2, IXHPOL2      Low & high screen coords for poloidal view x-axis
c  JYBPOL2, JYTPOL2      Bottom & top screen coords for poloidal view y-axis
c  IXLPLN2, IXHPLN2      Low & high screen coords for plan view x-axis
c  JYBPLN2, JYTPLN2      Bottom & top screen coords for plan view y-axis

c                        The following common block is for 3 plots per page!

      common /plot3/ ixlpol3, ixhpol3, jybpol3, jytpol3, ixlelv3, ixhelv3, 
     .   jybelv3, jytelv3, ixlpln3, ixhpln3, jybpln3, jytpln3

c  IXLPOL3, IXHPOL3      Low & high screen coords for poloidal view x-axis
c  JYBPOL3, JYTPOL3      Bottom & top screen coords for poloidal view y-axis
c  IXLELV3, IXHELV3      Low & high screen coords for elevation view x-axis
c  JYBELV3, JYTELV3      Bottom & top screen coords for elevation view y-axis
c  IXLPLN3, IXHPLN3      Low & high screen coords for plan view x-axis
c  JYBPLN3, JYTPLN3      Bottom & top screen coords for plan view y-axis

c                        The following common block is for the real coords of
c                          the plot box sizes

      common /plotr/ xlpol, xhpol, ybpol, ytpol, xlelv, xhelv, ybelv, ytelv, 
     .   xlpln, xhpln, ybpln, ytpln
      real*4 xlpol, xhpol, ybpol, ytpol, xlelv, xhelv, ybelv, ytelv, 
     .   xlpln, xhpln, ybpln, ytpln

c  XLPOL,   XHPOL        Low & high real coords (m) for poloidal view x-axis
c  YBPOL,   YTPOL        Bottom & top real coords (m) for poloidal view y-axis
c  XLELV,   XHELV        Low & high real coords (m) for elevation view x-axis
c  YBELV,   YTELV        Bottom & top real coords (m) for elevation view y-axis
c  XLPLN,   XHPLN        Low & high real coords (m) for plan view x-axis
c  YBPLN,   YTPLN        Bottom & top real coords (m) for plan view y-axis

c                        Text position in orbit plot

      common /plttxt/ ixtxt, jytxt
      integer ixtxt, jytxt

c  IXTXT                 Screen x-coord for left side of parameter listing
c  JYTXT                 Screen y-coord for top of parameter listing

c                        Color codes for plot elements follow:

      common /pltclr/ icfr, icorb, iclim, icflx, ictxt

c  ICFR                  Code for color of axes of graphs (defaults to black?)
c  ICORB                 Code for color of orbit
c  ICLIM                 Code for color of limiter
c  ICFLX                 Code for color of flux contours
c  ICTXT                 Code for color of text info

      common /psix/ npsiplt, simin, simax
      integer npsiplt
      real*8 simin, simax

c  NPSIPLT    number of psi contours to plot in output
c  SIMIN      minimum value of psi (poloidal flux) on grid
c  SIMAX      maximum value of psi (poloidal flux) on grid

      integer*4 dum
      integer iret

c     DUM         dummy variable to read comment lines in input file

c***************************begin executable code *******************

c FNAME is the parameter file name
      print *,  'rdpar: filename = ', FNAME 

c      OPEN(UNIT=23,FILE=FNAME,STATUS='old',ERR=1001, IOSTAT=iret)
      OPEN(UNIT=23,FILE=FNAME,STATUS='old')      

      PRINT 19,EFFLIM
19    FORMAT(' EFFLIM ',E10.2)


      READ(23,'(A80)') COMENT
d      write (8,*) 'rdpar comment = '// coment
      print *,  'rdpar: comment = ',coment


      READ(23, '(A4)') DUM
      print *,  'rdpar: DUM = ',DUM


c	  Now read name of file containing EFIT data
      READ(23, '(a80)') ifname
d      write (8,*) ' ifname=', ifname


c     Read name of file containing limiter description
      READ(23,'(a80)') limfnam
d      write (8, *) ' limfnam=', limfnam
      print *, 'rdpar : limfnam = ', limfnam

c     Read name of file containing limiter description
      READ(23,'(a80)') limfdir
d      write (8, *) ' limfdir=', limfdir
      print *, 'rdpar : limfdir = ', limfdir

c     Read name of file containing limiter description
      READ(23,'(a80)') ifdir
d      write (8, *) ' ifdir=', ifdir
      print *, 'rdpar : ifdir = ', ifdir 
     

c	  PORT is angle of port normal in R-Z plane
      READ(23, '(g16.7)') PORT
d      write (8,*) ' port=', port


      IF (eq99(PORT)) then
        ITERS(1)=1
        READ(23, '(g16.7)') PORT
        READ(23, '(g16.7)') STOP(1)
        READ(23, '(g16.7)') STEP(1)
      endif

111   continue


c	  PORTPH is angle of port normal wrt R-Z plane
      READ(23, '(g16.7)') portph


c     AL0 is initial gyroangle of orbit
      READ(23, '(g16.7)') AL0
d      write (8,*) ' portph=', portph, ' al0=', al0


c     if AL0 is -99.99 read the next 3 numbers as start, stop and step values for a range
c     of gyro angles
      IF (eq99(AL0)) then
        ITERS(3)=1
        READ(23, '(g16.7)') AL0
        READ(23, '(g16.7)') STOP(3)
        READ(23, '(g16.7)') STEP(3) 
      endif


c	  BE0 is initial "pitch angle" of orbit wrt toroidal direction
      READ(23, '(g16.7)') BE0
d      write (8,*) ' be0=', be0


      IF (eq99(BE0)) then
        ITERS(2)=1
        READ(23, '(g16.7)') BE0
        READ(23, '(g16.7)') STOP(2)
        READ(23, '(g16.7)') STEP(2)
      endif


c	  RD is major radius of detector, in cm.
c	  ZD is Z-position of detector, in cm.
c	  PHD is toroidal angle of detector, in radians.
c	  ZDET is charge of particle (in electron charges)
c	  ADET is mass of particle (in proton masses)
c	  ENER is energy of particle in MeV
      READ(23, '(g16.7)') RD
      READ(23, '(g16.7)') ZD
      READ(23, '(g16.7)') phd
      READ(23, '(g16.7)') ZDET,ADET
      READ(23, '(g16.7)') ENER
d      write (8,*) ' rd=',rd, ' zd=', zd, ' phd=', phd
d      write (8,*) ' zdet=', zdet, ' zdet=', zdet, ' ener=', ener


      IF (eq99(ENER)) then
        ITERS(5)=1
        READ(23, '(g16.7)') ENER
        READ(23, '(g16.7)') STOP(5)
        READ(23, '(g16.7)') STEP(5)
      endif


c	  S is single step size, in m.
c	  SSTP is maximum orbit length, in m.
c	  TOL is orbit integrator tolerance param, (1e-5)
c	  IS is profile exponent for computing efficiency
c	  D is spacing between collimating aperture
c	  and detector aperture, in cm.
      READ(23, '(g16.7)') S,SSTP
      READ(23, '(g14.7)') tol
d     write (8,*) ' s=',s, ' sstp=', sstp, ' tol=', tol


      READ(23, '(A4)') DUM
      READ(23, '(A4)') DUM
      READ(23, '(A4)') DUM
      READ(23, '(F8.2)') is
      READ(23, '(g16.7)') D
d     write (8,*) ' is=', is, ' d=', d


      if (d .le. 0) then
        write (6, *) ' ** Spacing between detector & collimator is <0: '
     .    , d
        close(23)
        stop
      endif


c	  XC is half-width of colimator, in cm.
c	  YC is half-height of collimator, in cm
c	  XD is half-width of detector aper, in cm
c	  YD is half-height of detector aper, in cm
      READ(23, '(g16.7)') XC
      READ(23, '(g16.7)') YC
      READ(23, '(g16.7)') XD
      READ(23, '(g16.7)') YD
d     write (8,*) ' xc, yc=', xc, yc, ' xd, yd=', xd, yd


      if ((xc+xd .eq. 0.0) .or. (yc+yd .eq. 0.0)) then
        write (6,*) ' ** Detector & collimator dimensions in X or Y ',
     .     'sum to zero!'
        close (23)
      endif


c	  NX & NY subdivide the aperture in pieces
c	  IFOR is 0 for time-reversed orbit, 1 o.w.
c	  IPL1 is 0 for no orbit plot, 1,2,3=# plots
c	  IRMN is 1 to plot minimum minor radius of 
c	  each orbit, 0 o.w.
      READ(23, '(I10)') NX,NY
      READ(23, '(A4)') DUM
      READ(23, '(A4)') DUM
      READ(23, '(A4)') DUM
      READ(23, '(I10)') IFOR
      READ(23, '(I10)') IPL1
      READ(23, '(I10)') IRMN
      READ(23, '(I10)') iflxplt
d       write (8,*) ' nx, ny=', nx, ny, ' ifor=', ifor, ' ipl1=', ipl1
d       write (8,*) ' irmn=', irmn, ' iflxplt=', iflxplt


      READ(23, '(A4)') DUM
      READ(23, '(A4)') DUM
      READ(23, '(A4)') DUM


c	  V0 is major radius of vessel center
c	  RMAX has something to do with DIII-D vacuum
c	  vessel shape (not used in other versions)
      READ(23, '(g16.7)') v0
      READ(23, '(g16.7)') rmax
d       write (8, *) ' v0=', v0, ' rmax=', rmax


c	  W1 & w2 are plotting window sizes, in m
      READ(23, '(g16.7)') W1,W2
d        write (8, *) ' w1, w2=', w1, w2


      READ(23, '(A4)') DUM
      READ(23, '(A4)') DUM
      READ(23, '(A4)') DUM


c	  Read info for setting up 2 plots per page
      READ(23, '(I10)') ixlpol2
      READ(23, '(I10)') ixhpol2
      READ(23, '(I10)') jybpol2
      READ(23, '(I10)') jytpol2
      READ(23, '(I10)') ixlpln2
      READ(23, '(I10)') ixhpln2
      READ(23, '(I10)') jybpln2
      READ(23, '(I10)') jytpln2
d       write (8, *) ' ixlpol2, ixhpol2=', ixlpol2, ixhpol2, 
d    .  ' jybpol2, jytpol2=', jybpol2, jytpol2
d       write (8, *) ' ixlpln2, ixhpln2=', ixlpln2, ixhpln2, 
d    .  ' jybpln2, jytpln2=', jybpln2, jytpln2


      if (ixlpol2 .ge. ixhpol2) print *, ' **2/page poloidal view lower',
     .  ' x limit (screen coords) >= upper limit: ', ixlpol2, ixhpol2
      if (jybpol2 .ge. jytpol2) print *, ' **2/page poloidal view lower',
     .  ' y limit (screen coords) >= upper limit: ', jybpol2, jytpol2
      if (ixlpln2 .ge. ixhpln2) print *, ' **2/page plan view lower',
     .  ' x limit (screen coords) >= upper limit: ', ixlpln2, ixhpln2
      if (jybpln2 .ge. jytpln2) print *, ' **2/page plan view lower',
     .  ' y limit (screen coords) >= upper limit: ', jybpln2, jytpln2


      READ(23, '(A4)') DUM
      READ(23, '(A4)') DUM
      READ(23, '(A4)') DUM


c	  Read info for 3 plots per page
      READ(23, '(I10)') ixlpol3
      READ(23, '(I10)') ixhpol3
      READ(23, '(I10)') jybpol3
      READ(23, '(I10)') jytpol3
      READ(23, '(I10)') ixlelv3
      READ(23, '(I10)') ixhelv3
      READ(23, '(I10)') jybelv3
      READ(23, '(I10)') jytelv3
      READ(23, '(I10)') ixlpln3
      READ(23, '(I10)') ixhpln3
      READ(23, '(I10)') jybpln3
      READ(23, '(I10)') jytpln3
d       write (8, *) ' ixlpol3, ixhpol3=', ixlpol3, ixhpol3, 
d    .  ' jybpol3, jytpol3=', jybpol3, jytpol3
d       write (8, *) ' ixlelv3, ixhelv3=', ixlelv3, ixhelv3, 
d    .  ' jybelv3, jytelv3=', jybelv3, jytelv3
d       write (8, *) ' ixlpln3, ixhpln3=', ixlpln3, ixhpln3, 
d    .  ' jybpln3, jytpln3=', jybpln3, jytpln3


      if (ixlpol3 .ge. ixhpol3) print *, ' **3/page poloidal view lower',
     .  ' x limit (screen coords) >= upper limit: ', ixlpol3, ixhpol3
      if (jybpol3 .ge. jytpol3) print *, ' **3/page poloidal view lower',
     .  ' y limit (screen coords) >= upper limit: ', jybpol3, jytpol3
      if (ixlelv3 .ge. ixhelv3) print *, ' **3/page elevation view lower',
     .  ' x limit (screen coords) >= upper limit: ', ixlelv3, ixhelv3
      if (jybelv3 .ge. jytelv3) print *, ' **3/page elevation view lower',
     .  ' y limit (screen coords) >= upper limit: ', jybelv3, jytelv3
      if (ixlpln3 .ge. ixhpln3) print *, ' **3/page plan view lower',
     .  ' x limit (screen coords) >= upper limit: ', ixlpln3, ixhpln3
      if (jybpln3 .ge. jytpln3) print *, ' **3/page plan view lower',
     .  ' y limit (screen coords) >= upper limit: ', jybpln3, jytpln3


      READ(23, '(A4)') DUM
      READ(23, '(A4)') DUM
      READ(23, '(A4)') DUM


c	  Read info for real extent of coordinates (m)
      READ(23, '(g16.7)') xlpol
      READ(23, '(g16.7)') xhpol
      READ(23, '(g16.7)') ybpol
      READ(23, '(g16.7)') ytpol
      READ(23, '(g16.7)') xlelv
      READ(23, '(g16.7)') xhelv
      READ(23, '(g16.7)') ybelv
      READ(23, '(g16.7)') ytelv
      READ(23, '(g16.7)') xlpln
      READ(23, '(g16.7)') xhpln
      READ(23, '(g16.7)') ybpln
      READ(23, '(g16.7)') ytpln
d       write (8, *) ' xlpol, xhpol=', xlpol, xhpol, 
d    .  ' ybpol, ytpol=', ybpol, ytpol
d       write (8, *) ' xlelv, xhelv=', xlelv, xhelv, 
d    .  ' ybelv, ytelv=', ybelv, ytelv
d       write (8, *) ' xlpln, xhpln=', xlpln, xhpln, 
d    .  ' ybpln, ytpln=', ybpln, ytpln


      if (xlpol .ge. xhpol) print *, ' **Poloidal view lower',
     .  ' x limit (real coords) >= upper limit: ', xlpol, xhpol
      if (ybpol .ge. ytpol) print *, ' **Poloidal view lower',
     .  ' y limit (real coords) >= upper limit: ', ybpol, ytpol
      if (xlelv .ge. xhelv) print *, ' **Elevation view lower',
     .  ' x limit (real coords) >= upper limit: ', xlelv, xhelv
      if (ybelv .ge. ytelv) print *, ' **Elevation view lower',
     .  ' y limit (real coords) >= upper limit: ', ybelv, ytelv
      if (xlpln .ge. xhpln) print *, ' **Plan view lower',
     .  ' x limit (real coords) >= upper limit: ', xlpln, xhpln
      if (ybpln .ge. ytpln) print *, ' **Plan view lower',
     .  ' y limit (real coords) >= upper limit: ', ybpln, ytpln


      READ(23, '(A4)') DUM
      READ(23, '(A4)') DUM
      READ(23, '(A4)') DUM


c	  Read info for colors in plots
      READ(23, '(I10)') ixtxt
      READ(23, '(I10)') jytxt
d       write (8, *) ' ixtxt=', ixtxt, ' jytxt=', jytxt


      READ(23, '(A4)') DUM
      READ(23, '(A4)') DUM
      READ(23, '(A4)') DUM


c	  Read info for colors in plots
      READ(23, '(I10)') icfr
      READ(23, '(I10)') icorb
      READ(23, '(I10)') iclim
      READ(23, '(I10)') icflx
      READ(23, '(I10)') ictxt
d       write (8, *) ' icfr=', icfr, ' icorb=', icorb, ' iclim=', iclim,
d    .  ' icflx=', icflx, ' ictxt=', ictxt


      READ(23, '(A4)') DUM
      READ(23, '(A4)') DUM
      READ(23, '(A4)') DUM


c	  Read number of flux contours to plot in poloidal view
      READ(23, '(I10)') npsiplt


      CLOSE(UNIT=23)
      return

c	  execute below if file open error

c 1001  CALL ERRSNS(I,J)
1001  PRINT 1011,IRET
1011  FORMAT(' PARAMETER INPUT FILE OPEN ERROR:',I6,2X/
     1  ' SEE COGNIZANT PROGRAMMER.')
      call finitt (0, 0)

      stop


      end

c
c --------------------------------- eq99 ------------------------------------
c return .true if the argument is -99.99
c     
      logical function eq99(d)
      real*8 d,e
      e = abs(d+99.99)
      eq99 = e.lt..0001d0
      end

