c
C     implicit integer*4 (i-n), real*8 (a-h, o-z)
      parameter (nfcoil=18,nsilop=41,nrogow=1,ntime=1,nacoil=1)
c-------------------------------------------------------------------
c--  New E-coil connections                       LLao, 95/07/11  --
c-------------------------------------------------------------------
      parameter (nesum=6)
      parameter (magpri67=29,magpri322=31,magprit=6)
      parameter (magpri=magpri67+magpri322)
      parameter (necoil=122,mpress=132,nvesel=24,ndata=41)
      parameter (nwwcur=18)
      parameter (nstark=16)
      parameter (nffcur=18,nppcur=18,npcurn=nffcur+nppcur
     .     ,mfnpcr=nfcoil+npcurn+nvesel+nwwcur+nesum+nfcoil
     .     ,npcur2=npcurn*2
     .     ,nrsmat=nsilop+magpri+nrogow+nffcur+1+npcurn+nwwcur+
     .      mpress+nfcoil+nstark
     .     ,nwcurn=nwwcur+npcurn,npcur3=npcurn*2
     .     ,nwcur2=nwcurn*2)
      parameter (npoint=800)
      parameter (nw=300,nh=300,nwnh=nw*nh)
      parameter (nh2=2*nh,nwrk=2*(nw+1)*nh)
      parameter (ncurrt=nvesel+nesum+nfcoil)
      parameter (mbdry=1500)
      parameter (nbwork=nsilop)
      parameter (kxiter=250,mqwant=30)
      parameter (nlimit=150,nlimbd=6)
      parameter (msbdry=mbdry+nsilop+nfcoil+1,msbdr2=2*msbdry)
      parameter (nrsma2=2*nrsmat)
      parameter (nwwf=2*nw)
      parameter (nwf=nwwf)
      parameter (nxtram=10,nxtrap=npoint)
      parameter (nxtlim=9,nco2v=3,nco2r=1)
      parameter (nangle=64)
      parameter (nfbcoil=12)
      parameter (kubicx = 4, kubicy = 4, lubicx = nw - kubicx + 1,
     .           lubicy = nh - kubicy + 1,
     .           kujunk = kubicx*kubicy*lubicx*lubicy)
      parameter (modef=4, modep=4, modew=4 , kubics=4 )
      parameter (nrsp=200)

      parameter (nxknot = lubicx+2*kubicx-1)
      parameter (nyknot = lubicy+2*kubicy-1)

