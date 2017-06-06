        SUBROUTINE MAGFLD
C***************************************************************
C***************************************************************
C   THIS ROUTINE CALCULATES THE MAGNETIC FIELD B AT R.  IT ALSO
C   CALCULATES THE  GYROFREQUENCY OMEGA.
C   IT USES THE FIELDS FROM THE SUBROUTINE PLOTMED.
C***************************************************************
C***************************************************************

        include 'dcalc.h'

        save

	
        common/consts/pm,echrg,emass,pi
        common/magcom/icur,cur,bt, scale
        common/geocom/rmax,v0,w1,w2
        common/p14com/zdet,ADET,ENER
        
        common/switch/ifor,ipl1,ipl2,irmn,rarray(201),rmn,ipoldir

C...B IS IN tesla. CUR IN KAMP.  ENER IN MEV.

        CALL PLOTMED(br,bz,bphi,absb)

        B(1)=br*scale
        B(3)=bz*scale
        B(2)=bphi*scale
        B(4)=absb*scale

        IF(IFOR.EQ.1) RETURN
		
c for ppro_1.INP ifor=0 so it should continue		

C...IF TIME-REVERSAL, THEN CHANGE SIGN OF MAGNETIC FIELD.

        DO 50 I=1,3
          B(I)=-B(I)
 50     CONTINUE

        RETURN
        END


