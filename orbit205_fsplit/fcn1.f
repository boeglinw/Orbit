


        SUBROUTINE FCN1(N,T,RV,RVPR)
C****************************************************************
C****************************************************************

      include 'dcalc.h'

      save

     

      DIMENSION RV(6),RVPR(6)
      real*8 bo(4)

c  BO  magnetic field components times omega

      common /fcncom/ omega

c	omega is the particle gyrofrequency divided by the local magnetic
c	field and the magnitude of the velocity.  It is, in fact, a 
c	constant of the motion, since it depends only on the particle's
c	energy and fundamental constants
		
        DO 10 I=1,3
           R(I)=RV(I)
10	continue
        R(4)=R(1)*r(1)+R(3)*r(3)

        CALL magfld

        DO 15 I=1,4
          BO(I)=B(I)*OMEGA
15      CONTINUE

        RVPR(1)=RV(4)
        RVPR(2)=RV(5)/RV(1)
        RVPR(3)=RV(6)
        RVPR(4)=RV(5)*rv(5)/RV(1)+RV(5)*BO(3)-RV(6)*BO(2)
        RVPR(5)=-RV(4)*RV(5)/RV(1)+RV(6)*BO(1)-RV(4)*BO(3)
        RVPR(6)=RV(4)*BO(2)-RV(5)*BO(1)

        RETURN
        END
