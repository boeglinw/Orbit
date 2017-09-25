      program read_eqdsk
c
      character*10 case_t(6) 
      parameter (MW=150)
      parameter (MH=150)
      parameter (MB_PTS = 1500)
      parameter (ML_PTS = 1500)
      parameter (neqdsk = 20)
      parameter (neqdsk_new = 26)

      dimension psirz(MW,MH),fpol(MW),pres(MW),ffprim(MW),
     >     pprime(MW),qpsi(MW),rbbbs(MB_PTS),zbbbs(MB_PTS), 
     >     rlim(ML_PTS),zlim(ML_PTS)

c
      
      open(neqdsk, file = 
     >     'eqdsk_old.dat', 
     >     status = 'old')

c      open(neqdsk_new, file = 
c     >     'eqdsk_new.dat', 
c     >     status = 'unknown')

      read (neqdsk,2000) (case_t(i),i=1,6),idum,nw,nh
D      print *, '(case_t(i),i=1,6),idum,nw,nh'
      write (*,2000) (case_t(i),i=1,6),idum,nw,nh

      read (neqdsk,2020) rdim,zdim,rcentr,rleft,zmid 
D      print *, 'rdim,zdim,rcentr,rleft,zmid'
      write (*,2020) rdim,zdim,rcentr,rleft,zmid 

c change sign
      read (neqdsk,2020) rmaxis,zmaxis,simag,sibry,bcentr 
D      print *, 'rmaxis,zmaxis,simag,sibry,bcentr'
      sibry = -1.*sibry
      simag = -1.*simag
      write (*,2020)rmaxis,zmaxis,simag,sibry,bcentr 
      
      read (neqdsk,2020) current,simag,xdum,rmaxis,xdum 
D      print *, 'current,simag,xdum,rmaxis,xdum'
c change sign
      simag = -1.*simag
      write (*,2020) current,simag,xdum,rmaxis,xdum 

      read (neqdsk,2020) zmaxis,xdum,sibry,xdum,xdum
c change sign
      sibry = -1.*sibry
D      print *, 'zmaxis,xdum,sibry,xdum,xdum'
      
      write (*,2020) zmaxis,xdum,sibry,xdum,xdum

      read (neqdsk,2020) (fpol(i),i=1,nw) 
D      print *,'(fpol(i),i=1,nw)'
      write (*,2020) (fpol(i),i=1,nw) 

      read (neqdsk,2020) (pres(i),i=1,nw) 
D      print *, '(pres(i),i=1,nw)'
      write (*,2020) (pres(i),i=1,nw) 

c change sign
      read (neqdsk,2020) (ffprim(i),i=1,nw) 
D      print *, '(ffprim(i),i=1,nw)'
D      print *, ' before : '
D      write(*,2020) (ffprim(i),i=1,nw)
      do i = 1, nw
         ffprim(i) = -1.*ffprim(i)
      enddo
D      print*, 'after :'
      write(*,2020) (ffprim(i),i=1,nw)


c change sign
      read (neqdsk,2020) (pprime(i),i=1,nw) 
D      print *, '(pprime(i),i=1,nw)'
D      print *, 'before :'
D      write(*,2020) (pprime(i),i=1,nw)
D      print*, 'after :'
      do i = 1, nw
         pprime(i) = -1.*pprime(i)
      enddo
      write(*,2020) (pprime(i),i=1,nw)


c change sign
      read (neqdsk,2020) ((psirz(i,j),i=1,nw),j=1,nh) 
D      print *, '((psirz(i,j),i=1,nw),j=1,nh)'
      do i = 1, nw
         do j = 1, nh
            psirz(i,j) = -1.*psirz(i,j)
         enddo
      enddo
      write (*,2020) ((psirz(i,j),i=1,nw),j=1,nh) 

c change sign
      read (neqdsk,2020) (qpsi(i),i=1,nw) 
D      print *, '(qpsi(i),i=1,nw) '
      do i = 1, nw
         qpsi(i) = -1.*qpsi(i)
      enddo
      write(*,2020) (qpsi(i),i=1,nw) 

      read (neqdsk,2022) nbbbs,limitr
D      print *,'nbbbs,limitr' 
      write(*,2022) nbbbs,limitr

      read (neqdsk,2020) (rbbbs(i),zbbbs(i),i=1,nbbbs) 
D      print *, '(rbbbs(i),zbbbs(i),i=1,nbbbs)'
      write(*,2020) (rbbbs(i),zbbbs(i),i=1,nbbbs)

      read (neqdsk,2020) (rlim(i),zlim(i),i=1,limitr)
D      print *, '(rlim(i),zlim(i),i=1,limitr)'
      write(*,2020) (rlim(i),zlim(i),i=1,limitr)
c 
2000  format (6a8,3i4) 
2020  format (5e16.9) 
2022  format (2i5)

      close (neqdsk)
      stop
      end
