       PROGRAM check_field
	   
C     **************************************************************
C read orbit input and chalculat the field for a given position
C     ************************************************************** 

      include 'dcalc.h'

      common/efitcom/ifname, ifdir
      common/limnam/limfnam, limfdir
      character*60 limfnam, limfdir

      common /iounit/iorbit, icounter, orbit_fname
      character *80 orbit_fname

      common/consts/pm,echrg,emass,pi
      common/magcom/icur,cur,bt
      common/geocom/rmax,V0,W1,W2
      common/p14com/zdet,adet,ener 
      parameter (mxorpt=64000)
      common/data5/x(mxorpt+1),y(mxorpt+1),z(mxorpt+1),
     .    x2(mxorpt+1),y2(mxorpt+1), b_r(mxorpt+1), b_phi(mxorpt+1), 
     .    b_z(mxorpt+1)
      common /data6/ iflxplt, steppsi(mxorpt+1), pphi(mxorpt+1),
     .    mu(mxorpt+1), energy(mxorpt+1), vpar(mxorpt+1), vperp(mxorpt+1),
     .    bmod(mxorpt+1), rho(mxorpt+1), lborho(mxorpt+1)
      real*8 mu, lborho

c  IFLXPLT    1 if particle flux surface vs step is to be plotted, 0 o.w.
c  STEPSPI    poloidal flux at particle position as a function of time
c  PPHI       canonical toroidal momentum of particle as function of time
c  MU         magnetic moment of particle as function of time
c  ENERGY     energy of particle as function of time
c  VPAR       parallel velocity of particle as function of time
c  VPERP      perpendicular velocity of particle as function of time
c  BMOD       magnitude of B at particle location as function of time
c  RHO        particle gyroradius as function of time
c  LBORHO     magnetic field scale length over particle gyroradius as fcn of t

c ORBIT_FNAME character string for the orbit data used for em. fitting
c IORBIT      unit number for orbit data

c  COMNT      character string from input file giving purpose of run
c  IFNAME     name of input files with EFIT magnetics data of plasma
c  LIMFNAM    name of file containing limiter configuration information
c  PORT       Orientation of the normal to the port in the R-Z plane
c  AL0        gyroangle in radians
c  BE0        complementary pitch angle in radians; +pi/2 is co-going
c  RD         major radius of detector
c  ZD         vertical position of detector
c  ZDET       charge of particle (in electron charges)
c  ADET       atomic mass of particle
c  ENER       energy of particle in MeV
c  S          length of a single step, in m. Does not affect accuracy of
c                integration
c  SSTP       maximum allowable length of orbit, in m.
c  IS         emission profile exponent (psi**is is emission function)
c  D          collimator-detector spacing, in m.
c  XC         collimator half-length, in m.
c  YC         collimator half-width, in m.
c  XD         detector half-length, in m.
c  YD         detector half-width, in m.
c  NX         number of divisions to use across length of detector
c  NY         number of divisions to use across width of detector
c  IFOR       0=time-reversed orbit computed; 1=forward orbit
c  IPL1       number of views of orbit to plot--0=no plot
c  IRMN       0=don't call RMIN, 1=do call RMIN
c  IFLXPLT    1 if particle flux surface vs step is to be plotted, 0 o.w.
c  STEPSPI    poloidal flux at particle position as a function of time
c  PPHI       canonical toroidal momentum of particle as function of time
c  MU         magnetic moment of particle as function of time
c  ENERGY     energy of particle as function of time
c  W1         plot window 1 width, in m.
c  W2         plot window 2 width, in m.
   
      common/shftcm/is
      common/switch/ifor,ipl1,ipl2,irmn,rarray(201),rmn,ipoldir
      common/char/fname,coment
      common /versn/ codenm, compdt, comptm

      character*4 name(5)
      integer iters(5)
      real*8 start(5)
      real*8 stop(5)
      real*8 step(5)

      common/iter/iters,start,stop,step,name

      dimension save(5) 
      real*8 is
      character*80 ifname 
      character*80 ifdir
      character*80 ifil
      character*80 fname

      real*8 rr, zz, p_phi
	  
c  coment    comment field describing this run of code
c  codenm    name of this version of code
c  compdt    compilation date of code (manual for now, auto later?)
c  comptm    compilation time of code (")

      character*80 coment
      character*9  compdt, comptm
      character*35 codenm

      logical found_one, outside, continue_iteration
      logical ifread, use_namelist

c  PM 'proton mass' in original version of code, but is actually
c  the mass (in kg) of 1 amu.  This and the other physical
c  constants below are from the NIST web site, 
c                    http://physics.nist.gov/cuu/Constants/
	  
      data pm/1.66053873d-27/, echrg/1.602176462d-19/,
     +  emass/9.10938188d-31/, pi/3.141592653590/

c  set up parameters for NSTX vacuum vesel and plasma

c      data (name(i),i=1,5) /'port','be0','al0','is','ener'/      
      data name(1) /'port'/
      data name(2) /'be0'/	  
      data name(3) /'a10'/
      data name(4) /'is'/
      data name(5) /'ener'/
      
      data codenm /'PPPL Lorentz ORBIT v205 for NSTX   '/
      data compdt /'00/07/28'/, comptm /'11:55:00'/

c orbit data file unit number       
      iorbit = 55

d     open(8, file = 'orbit205_dbg.dat', status = 'new')


C     INITIALIZE ITERATIVE CODES.

      DO I=1,5
         ITERS(I)=0
         START(I)=0.
	 STOP(I)=0. 
	 STEP(I)=0.
      enddo
      RARRAY(1)=0.0
      J=0
C initialize counter
      icounter = 0
C     READ IN PARAMTERS FROM DATA FILE.

      call chread('ENTER PARAMETER INPUT FILE NAME :', ifil, i1,i2)
      fname = ifil(i1:i2)
      use_namelist = .True.
      
      PRINT *, ' USING PARAMETER FILE: ', ifil, 
     >     ' using namelist = ',use_namelist 

C     READ INPUT PARAMETERS.

d     print *, "will call rdpar"

      call rdpar_nml

d     print *, "1.ITERS = ", ITERS

C     READ MAGNETICS DATA AND SET UP SPLINES.

      call pmedinit
d     print *, "2.ITERS = ", ITERS

C     READ FILE CONTAINING DESCRIPTION OF LIMITER.
      
      call liminit
d     print *, "3.ITERS = ", ITERS

      IPL2 = 0
      START(4)=IS
      START(5)=ENER
      
      
      do j = 1, detec
           port = ports(j)
           portph = portsph(j)
           RD = rds(j)
           be0 = be0s(j)
           al0 = al0s(j)
           zd = zds(j)
           phd = phda(j)
           call chrien(j)
           
        end do
        
!     d     print*,'called chrien in orbit205 iters(1) == 1'
        
        CLOSE(UNIT=24)
D     close (unit=8)

        rr = 0.
        
        do while (rr .ge. 0. )
           print *, ' Enter r,z,phi position : '
           read(*,*) rr, zz, pp
           r(1) = rr
           r(2) = pp
           r(3) = zz

           call magfld
           
           br = B(1)
           bz = B(3)
           bphi = B(2)
           bmag = B(4)
           
           print *, 'Field: Br, Bz, Bphi = ', br, bz, bphi
        enddo
c     close orbit data file      
        close(53)
        
        STOP 
        
        END 
