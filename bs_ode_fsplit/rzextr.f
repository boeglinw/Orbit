

c%
c======================================================/rzextr.for
c From Numerical Recipes, ch 15---Bulirsch-Stoer Method

c Subroutine to use a diagonal function extrapolation to evaluate NV functions
c at X=0 by fitting a diagonal rational function to a sequence of estimates
c with progressively smaller values X=XEST, and corresponding function vectors
c YEST.  This call is number IEST in the sequence of calls.  The extrapolation
c uses at most the last NUSE estimates.  Extrapolated function values are
c output as YZ, and their estimated error is output as DY.

c Rational function extrapolation can fail.  The extrapolated function might 
c have a pole at the desired evaluation point.  Or, more commonly, there might
c two poles that very nearly cancel, so that roundoff becomes a problem.  In
c this routine, the test for division by zero prevents program crash, but
c disguises the fact that the quantity computed as B1-C may have lost all
c significance.



      SUBROUTINE RZEXTR(IEST,XEST,YEST,YZ,DY,NV,NUSE)

      implicit none

      save
      integer*4 iest, nv, nuse
      real*8 xest, yest(nv), yz(nv), dy(nv)



c local declarations
      integer*4 imax, nmax, ncol, j, m1, k
      PARAMETER (IMAX=11,NMAX=10,NCOL=7)
      real*8 d(nmax,ncol), fx(ncol), c, yy, v, b1, b, ddy
      real*8 x(imax)


c------------------- executable --------------------



c save current independent value

      X(IEST)=XEST

      IF(IEST.EQ.1) THEN	!first try
        DO 11 J=1,NV
          YZ(J)=YEST(J)
          D(J,1)=YEST(J)
          DY(J)=YEST(J)
11      CONTINUE
      ELSE
        M1=MIN(IEST,NUSE)	!use at most NUSE previous members
        DO 12 K=1,M1-1
          FX(K+1)=X(IEST-K)/XEST
12      CONTINUE
        DO 14 J=1,NV		!evaluate next diagonal in tableau
          YY=YEST(J)
          V=D(J,1)
          C=YY
          D(J,1)=YY
          DO 13 K=2,M1
            B1=FX(K)*V
            B=B1-C
            IF(B.NE.0.) THEN
              B=(C-V)/B
              DDY=C*B
              C=B1*B
            ELSE		!avoid divide by 0
              DDY=V
            ENDIF
            if(k.ne.m1) V=D(J,K)
            D(J,K)=DDY
            YY=YY+DDY	
13        CONTINUE
          DY(J)=DDY	!returned values
          YZ(J)=YY
14      CONTINUE
      ENDIF
      RETURN
      END
