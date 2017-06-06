
c%
c=================================================/pzextr.for
c Polynomial extrapolation--less powerful than rational function extrapolation
c Use when HDID < HTRY (stepsize accomplished < desired stepsize)
c See Numerical Recipes, ch 15.

      SUBROUTINE PZEXTR(IEST,XEST,YEST,YZ,DY,NV,NUSE)

      implicit none
      save
      integer*4 iest, nv, nuse
      real*8 xest, yest(nv), yz(nv), dy(nv)

      integer*4 imax, ncol, nmax, j, m1, k1
      PARAMETER (IMAX=11,NCOL=7,NMAX=10)

      real*8 qcol(nmax,ncol), x(imax), d(nmax)
      real*8 delta, f1, f2, q


c------------------------ executable --------------------



      X(IEST)=XEST	!save current independent value

      DO 11 J=1,NV
        DY(J)=YEST(J)
        YZ(J)=YEST(J)
11    CONTINUE
      IF(IEST.EQ.1) THEN	!store first estimate in first column
        DO 12 J=1,NV
          QCOL(J,1)=YEST(J)
12      CONTINUE
      ELSE
        M1=MIN(IEST,NUSE)	!use at most nuse previous estimates
        DO 13 J=1,NV
          D(J)=YEST(J)
13      CONTINUE
        DO 15 K1=1,M1-1
          DELTA=1./(X(IEST-K1)-XEST)
          F1=XEST*DELTA
          F2=X(IEST-K1)*DELTA
          DO 14 J=1,NV		!propagate tableau 1 diagonal more
            Q=QCOL(J,K1)
            QCOL(J,K1)=DY(J)
            DELTA=D(J)-Q
            DY(J)=F1*DELTA
            D(J)=F2*DELTA
            YZ(J)=YZ(J)+DY(J)
14        CONTINUE
15      CONTINUE
        DO 16 J=1,NV
          QCOL(J,M1)=DY(J)
16      CONTINUE
      ENDIF
      RETURN
      END
