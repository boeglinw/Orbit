

        SUBROUTINE THAL4(N)
C***************************************************************
C***************************************************************
C   THIS VERSION TESTS WHETHER THE ION INTERSECTED THE VESSEL
C   WALLS. 
C***************************************************************
C***************************************************************
      include 'dcalc.h'

      save

      integer plotmed2

      
      COMMON/MAGCOM/ICUR,CUR,BT
      COMMON/GEOCOM/RMAX,V0

      common/consts/pm,echrg,emass,pi
      common/p14com/zdet,adet,ener 
      parameter (mxorpt=64000)
      common/switch/ifor,ipl1,ipl2,irmn,rarray(201),rmn,ipoldir
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

      real*8 mv0, hm, zre, ri0, pphisgn, vtotl2, vtmp2

c  MV0        mass of particle times its initial energy
c  HM         half the mass of the particle, in kg
c  ZRE        charge of particle in Coul
c  RI0        final number of steps computed in orbit (as real number)
c  PPHISGN    +1 or -1 depending upon whether orbit is computed forward or
c               backward in time (used in computation of P-phi)
c  VTOTL2     square of total velocity at each step in orbit
c  VTMP2      v-perp squared, but needs to be tested to see if >0

      common /fcncom/ omega
c
c	omega is the particle gyrofrequency divided by the local magnetic
c	field and the magnitude of the velocity.  It is, in fact, a 
c	constant of the motion, since it depends only on the particle's
c	energy and fundamental constants
c
      common /gridok/ bingrid
      logical bingrid

c     bingrid is set to .false. by PLOTMED if the particle is at a point
c     outside of the range of B-field data from EFIT.

      common /termm/ termmsg
      character*60 termmsg

c     termmsg contains a message generated in THAL4 to indicate why orbit
c     calculation was stopped.

      common /ptch/ riptch

c  RIPTCH	the real initial pitch angle relative to the magnetic field

      common /consv/ eneri, enerf, dener, mui, muf, dmu
      real*8 eneri, enerf, dener, mui, muf, dmu

c  ENERI, ENERF  initial and final energies of particle
c  DENER	     fractional change in energy of particle over orbit (ideally 0)
c  MUI, MUF  	 initial and final magnetic moments of particle
c  DMU		     fractional change in magnetic moment over orbit

      DIMENSION RV(6),C(24),W(6,9)
      EXTERNAL FCN1
      logical hitlim  

C     HITLIM is TRUE if particle hits limiter, FALSE otherwise.
	   
      bingrid = .true.
      RMAX2=RMAX**2.
	  
C   COMPUTE OMEGA ONCE, SINCE IT IS CONSTANT

      omega = zdet * echrg / sqrt(2.0d6 * echrg * ener * adet * pm)

      mv0 = pm * adet * sqrt(2.0d6 * echrg * ener / (adet * pm)) 
      hm = ener * echrg * 1.0d6
      zre = zdet * echrg 
      pphisgn = -1.0d0
      if (ifor .ge. 1) pphisgn = 1.0d0
		
c       ifor=0 for ppro_1.INPU so this if statement should not pass		

D       write (8, *) ' mv0=', mv0, ' hm=', hm, ' zre=', zre

      X(2)=R(1)
      Y(2)=R(2)
      Z(2)=R(3)
      X2(2)=R(1)*COS(R(2))
      Y2(2)=R(1)*SIN(R(2))
      NW=6
      N1=6
      T=0.0d0
      DO 10 I=1,3
         RV(I)=R(I)
         RV(I+3)=V(I)/V(4)
 10   continue
      V(4)=1.0d0
      IND=1
      DELT=S/V(4)

c     If flux surface postion plot is turned on, then compute ploidal
c     flux and conserved at start position of particle. 
      if (plotmed2 (flux, rv(1), rv(3)) .eq. 0) then
c        outside the boundary set flux to a huge number          
         flux = 1e10
      endif
      steppsi(2) = flux
      call magfld


        
c     compute the real initial pitch angle wrt B-field at starting point.
c     That  angle is the arccos of v dot B /|v||B|.

      if (b(4) .ne. 0.0d0) then
         vpar(2) = (v(1) * b(1) + v(2) * b(2) + v(3) * b(3)) / b(4)
         rpt = vpar(2)
         riptch = 1.80d2/pi*acos(rpt)
      else
         riptch = -99.99d0
      endif

D     write (8,*) ' rpt=', rpt, ' riptch=', riptch

      if (iflxplt .gt. 0) then
         bmod(2) = b(4)
         
         
C     NEED TO DEFINE CANONICAL MOMENTUM WITH NEGATIVE SIGN DUE TO
C     EFIT'S DEFITION OF THE FLUX, IF GOING BACKWARS IN TIME, ELSE
C     POSITIVE SIGN IF GOING FORWARDS IN TIME
         
         pphi(2) = mv0 * rv(1) * rv(5) + pphisgn * zre * flux 
         energy(2) = ener * v(4) * v(4)
         if (b(4) .ne. 0.0d0) then 
            vtmp2 = v(4) * v(4) - vpar(2) * vpar(2)
            if (vtmp2 .ge. 0.0d0) then
               vperp(2) = sqrt(vtmp2)
            else
               vperp(2) = 0.0d0
            endif
            mu(2) = hm * vtmp2 / b(4)
            rho(2) = vperp(2) / (omega * b(4))
            lborho(2) = 0.0d0
         else
            vpar(2) = 1.0d0
            vperp(2) = 0.0d0
            mu(2) = 0.0d0
            rho(2) = 0.0d0
            lborho(2) = -1.0d0
         endif
         mui = mu(2)
         
D     write (8, *) ' For initial position, flux=', flux, ' P-phi=',
D     .       pphi(2), ' energy=', energy(2), ' mu=', mu(2), 
D     .       ' vpar=', vpar(2), ' vperp=', vperp(2),
D     .       ' |B|=', b(4) , ' vtmp2=', vtmp2
         
      endif
      
c     write out position & velocity at each step
      
c     RV(1) : R
c     RV(2) : phi
c     RV(3) : Z          
c     B(1)  : br
c     B(2)  : bphi
c     B(3)  : bz
c     save initial magnetic field
      b_r(2) = b(1)
      b_phi(2) = b(2)
      b_z(2) = b(3)
      
      write (53, 4010) 1, (rv(jj), jj= 1, 6), (b(jj), jj = 1,3)
      
C     COMPUTE ALL THE STEPS IN THE ORBIT.
      
      DO 100 I=2,N
         I0=I
         TEND=REAL(I-1)*DELT
         ier = 0
         
C     12AUG97 DSD SEEKING TO REPLACE OLD IMSL ROUTINE DVERK WITH
C     BS_ODE, BORROWED FROM TFTR VERSION OF THIS CODE, AND
C     ORIGINATING IN NUMERICAL RECIPES.
         
         call bs_ode (fcn1, n1, rv, t, delt, tol, ier)
         if (ier .ne. 0) print *, ' BS_ODE returned error code ', ier
c     RV(1) : R
c     RV(2) : phi
c     RV(3) : Z          
c     B(1)  : br
c     B(2)  : bphi
c     B(3)  : bz
         X(I+1)=RV(1)
         Y(I+1)=RV(2)
         Z(I+1)=RV(3)
         X2(I+1)=RV(1)*COS(RV(2))
         Y2(I+1)=RV(1)*SIN(RV(2))
         
         b_r(I+1) = b(1)
         b_phi(I+1) = b(2)
         b_z(I+1) = b(3)
         
         RR=(RV(1)-V0) * (rv(1)-v0) + RV(3) * rv(3)
         
c     write out position & velocity and magnetic field at each step
c     RV(1) : R
c     RV(2) : phi
c     RV(3) : Z          
c     B(1)  : br
c     B(2)  : bphi
c     B(3)  : bz
         write (53, 4010) I, (rv(jj), jj= 1, 6), (b(jj), jj = 1,3)
         
C     CHECK TO SEE WHETHER PARTICLE ENTERED REGION FOR WHICH THERE IS
C     NO B-FIELD DATA.
         
         if (.not. bingrid) go to 103
         
C     IF FLUX SURFACE POSITION PLOT IS TURNED ON, THEN COMPUTER POLOIDAL
C     FLUX AT PRESENT POSITION OF PARTICLE. ALSO COMPUTER CONSERVED
C     QUANTITIES AND ACCUMULATE IN ARRAYS
         
         
         if (plotmed2 (flux, rv(1), rv(3)) .eq. 0) then
c     outside the boundary set flux to a huge number          
            flux = 1.e10
         endif
         steppsi(i+1) = flux
         
         if (iflxplt .gt. 0) then
            
C     NEED TO DEFINE CANONICAL MOMENTUM WITH NEGATIVE SIGN DUE TO EFIT'S
C     DEFINITION OF THE FLUX, IF GOING BACKWARDS IN TIME, ELSE POSITIVE
C     SIGN IF GOING FORWARDS IN TIME. 
            
            pphi(i+1) = mv0 * rv(1) * rv(5) + pphisgn * zre * flux 
            vtotl2 = rv(4) * rv(4) + rv(5) * rv(5) + rv(6) * rv(6)
            energy(i+1) = ener * vtotl2
            call magfld
            bmod(i+1) = b(4)
            if (b(4) .ne. 0.0d0) then 
               vpar(i+1) = (rv(4) * b(1) + rv(5) * b(2) +
     .              rv(6) * b(3)) / b(4)
               vtmp2 = vtotl2 - vpar(i+1) * vpar(i+1)
               if (vtmp2 .ge. 0.0d0) then
                  vperp(i+1) = sqrt(vtmp2)
               else
                  vperp(i+1) = 0.0d0
               endif
               mu(i+1) = hm * vperp(i+1) * vperp(i+1) / b(4)
            else
               vpar(i+1) = 1.0d0
               vperp(i+1) = 0.0d0
               mu(i+1) = 0.0d0
            endif
            rho(i+1) = vperp(i+1) / (omega * b(4))
            lborho(i+1) = 0.0d0
            
            
D     write (8, *) ' For step ', i, ', flux=', flux, ' P-phi=',
D     .          pphi(i+1), ' energy=', energy(i+1), ' mu=', mu(i+1), 
D     .          ' vpar=', vpar(i+1), ' vperp=', vperp(i+1),
D     .          ' |B|=', b(4), ' vtotl2=', vtotl2, ' vtmp2=', vtmp2,
D     .          ' rho=', rho(i+1)
            
         endif
         
C     CHECK FOR PARTICLE INTERSECTING LIMITER.
C     IF INTERSECT LIMITER, THEN PRINT MESSAGE AND END THE CALCULATION.
C     FILE DOES NOT EXIST FOR THIS SUBROUTINE.
         
         call intlim (hitlim)
         if (hitlim) go to 101
         
 100  continue
      
C     "NORMAL" TERMINATION HERE- HIT STEP LIMIT. SET MESSAGE TO PRINT.
      
      write (termmsg, 99) i0
 99   FORMAT('MAXIMUM NUMBER OF STEPS HAS BEEN REACHED, N=',I6,10x)		
      GOTO 110 
      
 101  continue
      
      
C     TERMINATE HERE IF PARTICLE HIT LIMITER.
      write (termmsg, 95) r(1), r(2)*1.80d2/pi, r(3)
 95   FORMAT('ORBIT STRUCK LIMITER AT R=',F6.3,' Ph=', f11.1,
     .     ' Z=',F7.3,3x)        
      go to 110
      
 103  continue
      
C     TERMINATION DUE TO NO B-FIELD DATA.
      
      write (termmsg, 96) r(1), r(2)*1.80d2/pi, r(3)
 96   format('ORBIT ENTERED NO B DATA REGION,R=',F6.3,' Ph=',f6.1,
     .     ' Z=',F7.3)
      
 110  continue
      
 4010 format (1x, i6, 9(2x, g14.7))
      
C     PRINT REASON FOR END OF ORBIT COMPUTATION, UPDATE NUMBER
C     OF STEPS, THEN QUIT.
      print *,termmsg
      ri0 = dble (i0)
      N=I0
      X(1) = ri0
      Y(1) = ri0
      Z(1) = ri0
      X2(1) = ri0
      Y2(1) = ri0
      b_r(1) = ri0
      b_phi(1) = ri0
      b_z(1) = ri0
      steppsi(1) = ri0
      pphi(1) = ri0
      mu(1) = ri0
      energy(1) = ri0
      vpar(1) = ri0
      vperp(1) = ri0
      bmod(1) = ri0
      rho(1) = ri0
      lborho(1) = ri0
      RETURN
      
      
      END
      
