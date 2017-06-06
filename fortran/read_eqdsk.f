      program read_eqdsk
c
      character*10 case_t(6) 
      parameter (MW=100)
      parameter (MH=100)
      parameter (MB_PTS = 1500)
      parameter (ML_PTS = 1500)
      parameter (neqdsk = 20)

      dimension psirz(MW,MH),fpol(MW),pres(MW),ffprim(MW),
     >     pprime(MW),qpsi(MW),rbbbs(MB_PTS),zbbbs(MB_PTS), 
     >     rlim(ML_PTS),zlim(ML_PTS)

      common /results/ case_t, nw, nh, rdim, zdmin,rcentr, rleft, zmid,
     >     rmaxis, zmaxis, simag, sibry, bcenter,
     >     current, sigmag,
     >     psirz,fpol,pres,ffprim,
     >     pprime,qpsi,rbbbs,zbbbs, 
     >     rlim,zlim, nbbbs,limitr
c
      open(neqdsk, file = 
     >     'g135445.00805.EFIT01.mds.corrected.qscale_1.00000', 
     >     status = 'old')

      read (neqdsk,2000) (case_t(i),i=1,6),idum,nw,nh
      print *, '(case_t(i),i=1,6),idum,nw,nh'
      write (*,*) (case_t(i),i=1,6),idum,nw,nh

      read (neqdsk,2020) rdim,zdim,rcentr,rleft,zmid 
      print *, 'rdim,zdim,rcentr,rleft,zmid'
      write (*,*) rdim,zdim,rcentr,rleft,zmid 

      read (neqdsk,2020) rmaxis,zmaxis,simag,sibry,bcentr 
      print *, 'rmaxis,zmaxis,simag,sibry,bcentr'
      write (*,*)rmaxis,zmaxis,simag,sibry,bcentr 
      
      read (neqdsk,2020) current,simag,xdum,rmaxis,xdum 
      print *, 'current,simag,xdum,rmaxis,xdum'
      write (*,*) current,simag,xdum,rmaxis,xdum 

      read (neqdsk,2020) zmaxis,xdum,sibry,xdum,xdum
      print *, 'zmaxis,xdum,sibry,xdum,xdum'
      write (*,*) zmaxis,xdum,sibry,xdum,xdum

      read (neqdsk,2020) (fpol(i),i=1,nw) 
      print *,'(fpol(i),i=1,nw)'
      write (*,'(1pe10.3)') (fpol(i),i=1,nw) 

      read (neqdsk,2020) (pres(i),i=1,nw) 
      print *, '(pres(i),i=1,nw)'
      write (*,'(1pe10.3)') (pres(i),i=1,nw) 

      read (neqdsk,2020) (ffprim(i),i=1,nw) 
      print *, '(ffprim(i),i=1,nw)'
      write(*,'(1pe10.3)') (ffprim(i),i=1,nw)

      read (neqdsk,2020) (pprime(i),i=1,nw) 
      print *, '(pprime(i),i=1,nw)'
      write(*,'(1pe10.3)') (pprime(i),i=1,nw)

      read (neqdsk,2020) ((psirz(i,j),i=1,nw),j=1,nh) 

      read (neqdsk,2020) (qpsi(i),i=1,nw) 
      print *, '(qpsi(i),i=1,nw) '
      write(*,'(1pe10.3)') (qpsi(i),i=1,nw) 

      read (neqdsk,2022) nbbbs,limitr
      print *,'nbbbs,limitr' 
      write(*,*) nbbbs,limitr

      read (neqdsk,2020) (rbbbs(i),zbbbs(i),i=1,nbbbs) 
      print *, '(rbbbs(i),zbbbs(i),i=1,nbbbs)'
      write(*,'(1pe10.3, 2x, 1pe10.3)') (rbbbs(i),zbbbs(i),i=1,nbbbs)

      read (neqdsk,2020) (rlim(i),zlim(i),i=1,limitr)
      print *, '(rlim(i),zlim(i),i=1,limitr)'
      write(*,'(1pe10.3, 2x, 1pe10.3)') (rlim(i),zlim(i),i=1,limitr)
c 
2000  format (6a8,3i4) 
2020  format (5e16.9) 
2022  format (2i5)

      close (neqdsk)
      stop
      end
