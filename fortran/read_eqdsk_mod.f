c
c can be imported into Python
C 
      subroutine read_eqdsk(eq_file)
c
      character*10 case_t(6) 
      character*80 eq_file
      parameter (MW=129)
      parameter (MH=129)
      parameter (MB_PTS = 1500)
      parameter (ML_PTS = 1500)
      parameter (neqdsk = 20)

      dimension psirz(MW,MH),fpol(MW),pres(MW),ffprim(MW),
     >     pprime(MW),qpsi(MW),rbbbs(MB_PTS),zbbbs(MB_PTS), 
     >     rlim(ML_PTS),zlim(ML_PTS)

      common /results/ case_t, nw, nh, rdim, zdim, rcentr, rleft, zmid,
     >     rmaxis, zmaxis, simag, sibry, bcentr,
     >     current, sigmag,
     >     psirz,fpol,pres,ffprim,
     >     pprime,qpsi,rbbbs,zbbbs, 
     >     rlim,zlim, nbbbs,limitr
c
      eq_file = adjustl(eq_file)
      open(neqdsk, file = eq_file(1:len_trim(eq_file)),
     >     status = 'old')

      read (neqdsk,2000) (case_t(i),i=1,6),idum,nw,nh
      read (neqdsk,2020) rdim,zdim,rcentr,rleft,zmid 
      read (neqdsk,2020) rmaxis,zmaxis,simag,sibry,bcentr 
      read (neqdsk,2020) current,simag,xdum,rmaxis,xdum 
      read (neqdsk,2020) zmaxis,xdum,sibry,xdum,xdum
      read (neqdsk,2020) (fpol(i),i=1,nw) 
      read (neqdsk,2020) (pres(i),i=1,nw) 
      read (neqdsk,2020) (ffprim(i),i=1,nw) 
      read (neqdsk,2020) (pprime(i),i=1,nw) 
      read (neqdsk,2020) ((psirz(i,j),i=1,nw),j=1,nh) 
      read (neqdsk,2020) (qpsi(i),i=1,nw) 
      read (neqdsk,2022) nbbbs,limitr
      read (neqdsk,2020) (rbbbs(i),zbbbs(i),i=1,nbbbs) 
      read (neqdsk,2020) (rlim(i),zlim(i),i=1,limitr)
c 
2000  format (6a8,3i4) 
2020  format (5e16.9) 
2022  format (2i5)

      close (neqdsk)
      return
      end
