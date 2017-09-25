      subroutine rdeqdsk(ifnamein,imfit,mw,mh,ier)

c**********************************************************************
c**                                                                  **
c**     MAIN PROGRAM:  MHD FITTING CODE                              **
c**                                                                  **
c**                                                                  **
c**     SUBPROGRAM DESCRIPTION:                                      **
c**          weqdsk reads  out the GAQ type eqdsk.                   **
c**                                                                  **
c**     CALLING ARGUMENTS:                                           **
c**                                                                  **
c**     REFERENCES:                                                  **
c**          (1)                                                     **
c**          (2)                                                     **
c**                                                                  **
c**     RECORD OF MODIFICATION:                                      **
c**          29/06/83..........first created                         **
c**          02/02/99..........modified by DSD for PPPL orbit code   **
c**                                                                  **
c**                                                                  **
c**   WB 2017 removed all ioption < 1 part of the code.              **
c**********************************************************************

      implicit real*8 (a-h, o-z)

      save
      include 'eparmdu129.h'

      parameter (npitch=8*nlimit)

      common/cwork3/c(kubicx,lubicx,kubicy,lubicy),wk(nwrk),
     .         copy(nw,nh),bkx(lubicx+1),bky(lubicy+1),
     .         lkx,lky
      common/mcom/rmaxis,zmaxis
      common/limmm/xlmin,xlmax,ylmin,ylmax,radum
      common/cpitch/ipitch,tanbn(npitch),angbn(npitch),sisibn(npitch)
      common/input1/icondn,itek,kdata,ico2,itrace,ierchk,iconvr,ixray
      common/cwork1/psirz(nw,nh)
      common/parame/volp(nw),pprime(nw),pres(nw),ffprim(nw),fpol(nw)
      common/consta/pi,tmu,twopi,ioption,ipress,idolim,isurf
      common/nio/nin,nout,ntty,nrsppc,nrspfc,nttyo,neqdsk,nffile,nsave
      common/limite/limitr,xlim(nlimit),ylim(nlimit)
      common/cffol/qpsi(nw),bfpol(nw),cfpol(nw),dfpol(nw),mwfpol,
     .             rbdry(mbdry),zbdry(mbdry),nbdry,xxxsi(nw),sifm,sifb
     .            ,rzero,bzero,piii,qout95,pcurrt(nwnh)
      common/ccase/case(6)
      common/vtor/presw(nw),preswp(nw),kvtor,rvtor,cwrmid(nw)
     .            ,prw(nw),wpsi(nwnh),presst(nwnh),pressw(nwnh)
      common/cerror/delerr
      common/gtable/rgrid(nw),zgrid(nh)
      common/cgrid/darea,drgrid,dzgrid,qmaxis,cratio,dfsqe
      common/cpsi/psi(nwnh),psibry,simag,sidif,xpsi(nwnh),eouter
      common /dsdmhd/ pr(nw), prp(nw), ffp(nw)

      common/switch/ifor,ipl1,ipl2,irmn,rarray(201),rmn,ipoldir
      
      character ifnamein*160,ifname*80,case*10,efitd65*32
      character nshot5*5,nshot6*6,nntime*4,vernum*6,nine*1
      dimension workk(nw),coils(nsilop),expmp2(magpri)
     .  ,xsi(nw),bfp(nw),cfp(nw),dfp(nw),pds(6),brsp(nrsp),prexp(1)
     .  ,bpp(nw),cpp(nw),dpp(nw),fwtsi(nsilop)
     .  ,bpr(nw),cpr(nw),dpr(nw),bqpsi(nw),cqpsi(nw),dqpsi(nw)

      dimension bprwp(nw),cprwp(nw),dprwp(nw)
      dimension bprw(nw),cprw(nw),dprw(nw),pressu(1)
      dimension curmid(nw)

c  CURMID Work array (originally argument of this routine)

      namelist/out1/ishot,itime,betap0,rzero,qenp,enp,emp,plasma,
     .     expmp2,coils,prexp,btor,rcentr,brsp,icurrt,rbdry,zbdry,
     .     nbdry,fwtsi,fwtcur,mxiter,nxiter,limitr,xlim,ylim,error,
     .     iconvr,ibunmn
      data efitd65/'PHYS_DATA:[D3PHYS.DIIID.EFITD65]'/

c  Stuff below is spliced in from old WEQDSK routine obtained from Bill Heidbrink

      data nin/11/,nout/10/,ntty/5/,nrsppe/25/,nrspfe/26/,meontr/35/
     .     ,lfile/36/,neqdsk/38/,nffile/40/,nsave/49/,kltype/-1/
     .     ,nttyo/6/
      data tmu/2.0e-07/,reentr/1.6955/

c-----------------------------Executable code begins here-----------------

d     print*,' Entered routine WEQDSK which is called from pmedinit'
d	  print*,' which is in plotmed205.f'
      pi = 3.1415926535897932d0
      ier = 0
      ixray=0
      ierchk=0
      itek=4
      twopi=2.0*pi
      ioption = 1

c called from pmedinit then ifnamein= ifname2= g + ifname= g+ well.= gwell.
c this ifname is stored from scanned input file in subroutine rdpar

      print *, '--> EQ File unit, name : ', neqdsk, ifnamein
      open(unit=neqdsk,file=ifnamein,status='old',err=2030)
      
      READ(neqdsk, '(6a8,3i4)') (case(i),i=1,6),imfit,mw,mh
d     write (8, 400) (case(i), i=1,6)
d 400  format(' Case=', 6a10)
d     print*,'case(1)',case(1)
d     print*,'case(2)',case(2)
d     print*,'case(3)',case(3)
d     print*,'case(4)',case(4)
d     print*,'case(5)',case(5)
d     print*,'case(6)',case(6)
d     print*,'imfit= ',imfit
d     print*,'mw= ',mw
d     print*,'mh= ',mh
      
      READ(neqdsk, '(5e16.9)') xdim,zdim,rzero,rgrid(1),zmid
      print*,'xdim= ',xdim
      print*,'ydim= ',ydim
      print*,'rzero= ',rzero 
      print*,'rgrid(1)= ',rgrid(1)
      print*,'zmid= ',zmid
      print*,'mw= ',mw
      print*,'mh= ',mh
      
      READ(neqdsk, '(5e16.9)') rmaxis,zmaxis,ssimag,ssibry,bzero

      
      if (ipoldir .gt. 0) then
         ssimag = -ssimag
         ssibry = -ssibry
      endif

d     print*,'rmaxis= ',rmaxis 
d     print*,'zmaxis= ',zmaxis 
d     print*,'mw= ',mw 
d     print*,'mh= ',mh
      
      READ(neqdsk, '(5e16.9)') piii,xdum,xdum,rmaxis,xdum
d     print*,'xdum= ',xdum
d     print*,'mw= ',mw
d     print*,'mh= ',mh
      
      READ(neqdsk, '(5e16.9)') zmaxis,xdum,sdum,xdum,xdum
d     print*,'zmaxis= ',zmaxis
d     print*,'mw= ',mw
d     print*,'mh= ',mh
      
      READ(neqdsk, '(5e16.9)') (fpol(i),i=1,mw)
d     print*,'fpol(1)= ',fpol(1), ' fpol(',mw,') = ', fpol(mw) 
d     print*,'mw= ',mw
d     print*,'mh= ',mh
      
      READ(neqdsk, '(5e16.9)') (pres(i),i=1,mw) 
d     print*,'PRES(1)= ',pres(1),' pres(',mw,') = ', pres(mw)  
d     print*,'mw= ',mw
d     print*,'mh= ',mh
      
c thsi would be ffprim stored in temp. array workk      
      READ(neqdsk, '(5e16.9)') (workk(i),i=1,mw) 
d     print*,'workk(1)= ',workk(1), ' workk(',mw,') = ', workk(mw)
d     print*,'mw= ',mw
d     print*,'mh= ',mh
      
      drgrid=xdim/float(mw-1)
      dzgrid=zdim/float(mh-1)
      do 200 i=1,mw
         rgrid(i)=rgrid(1)+(i-1)*drgrid
 200  continue
      do 220 i=1,mh
         zgrid(i)=zmid-zdim/2.+(i-1)*dzgrid
 220  continue
      darea=drgrid*dzgrid
      do 310 i=1,mw
c        twopi*tmu = mu0 in SI units      
         if (imfit.ge.0) ffprim(i)=-workk(i)/(twopi*tmu)
         if (imfit.lt.0) ffprim(i)=-workk(i)
 310  continue
c pprime data stored in temp. work array
      READ(neqdsk, '(5e16.9)') (workk(i),i=1,mw)
      do 315 i=1,mw
         pprime(i)=-workk(i)
 315  continue
      
d     print*,'Finished reading PPRIME'  
d     print*,'mw= ',mw
d     print*,'mh= ',mh

c      do i = 1, mw
c         READ(neqdsk, '(5e16.9)') ((psirz(i,j),j=1,mh)
c      enddo
      READ(neqdsk, '(5e16.9)') ((psirz(i,j),i=1,mw),j=1,mh)
c d     do j = 1, mh
c d        do i = 1, mw
c d           print *, i, j , psirz(i,j)
c d        enddo      
c d     enddo      
d     print*, 'psirz(1,1) = ', psirz(1,1), ' psirz(',mw,mh,') = ', psirz(mw,mh)
c
c     switch polodial field direction
c
      if (ipoldir .gt. 0) then
         do i=1,mw
            do j=1,mh
               psirz(i,j) = -psirz(i,j)
            enddo
         enddo
      endif
      
      do i=1,mw
         do j=1,mh
            kk=(i-1)*mh+j
            psi(kk)=psirz(i,j)
         enddo
      enddo
      simag=ssimag
      psibry=ssibry
      sifm=simag
      sifb=psibry
      

      ioption=1
d       print*,'Starting to read QPSI' 
d       print*,'mw= ',mw
         
        READ(neqdsk, '(5e16.9)') (qpsi(i),i=1,mw)
d       print*, 'qpsi(1) =  ', qpsi(1), ' qpsi(',mw,') = ', qpsi(mw)        

d       print*,'Starting to Read Limiter and Boundary'        
        READ(neqdsk, '(2i5)') nbdry,limitr
         
        READ(neqdsk, '(5e16.9)') (rbdry(i),zbdry(i),i=1,nbdry)
        READ(neqdsk, '(5e16.9)') (xlim(i),ylim(i),i=1,limitr)
c     READ(neqdsk,out1)

        xlmin=xlim(1)
        xlmax=xlmin
        ylmin=ylim(1)
        ylmax=ylmin
        do 22140 i=2,limitr
            xlmin=min(xlmin,xlim(i))
            xlmax=max(xlmax,xlim(i))
            ylmin=min(ylmin,ylim(i))
            ylmax=max(ylmax,ylim(i))
22140   continue
      close(unit=neqdsk)
      
      dxsi=1./float(mw-1)
      do 750 i=1,mw
         xsi(i)=(i-1)*dxsi
         xxxsi(i)=xsi(i)
 750  continue
      
c************************************
      
30900 close(unit=neqdsk)
      
      xguess=(rgrid(1)+rgrid(mw))/2.
      radum=(xguess+xlmin)/2.

c      piii is plasma current
      if (piii.le.-1.e3) then
         negcur=1
d     print*,'negcur= ',negcur
      else
         negcur=0
d     print*,'negcur= ',negcur
      endif
      

c------------------------------------------------------------------------------
c------------------------------------------------------------------------------
         
c      rangle=(xlmin+xlmax)/2.
      rangle=rcentr
c store psi back in 2d array copy
      do 1000 i=1,mw
         do 1000 j=1,mh
            kk=(i-1)*mh+j
            copy(i,j)=psi(kk)
 1000 continue
c----------------------------------------------------------------------
c--   fit 2-d zpline to psi                                          --
c----------------------------------------------------------------------
c     call ibcccu(copy,rgrid,mw,zgrid,mh,c,nw,wk,ier)
c     calculate a set of 2d-spline coefficients, psi is a 2d array mapped on 1d
c used later to interpolate on the grid
      call sets2d(psi,c,rgrid,mw,bkx,lkx,zgrid,mh,bky,lky,wk,ier)
      mwfpol=mw
c     calc a set of 1d splines coefficients for fpol as a function of a coordinate
c     where 0 is at the location of the magnetic axis and 1 at the last closed
c     flux surface: xsi
      call zpline(mwfpol,xxxsi,fpol,bfpol,cfpol,dfpol)
      if (ipitch.le.0) return
c
c calculate values of psi on limiter
      do 1005 i=2,limitr
         delx=(xlim(i)-xlim(i-1))/4.
         dely=(ylim(i)-ylim(i-1))/4.
         do 1005 k=1,4
            xwant=xlim(i-1)+(k-1)*delx
            ywant=ylim(i-1)+(k-1)*dely
c     call dbcevl(rgrid,mw,zgrid,mh,c,nw,xwant,ywant,pds,ier) 
            call seva2d(bkx,lkx,bky,lky,c,mw,mh,xwant,ywant,pds,ier,n333)
            abpolz=pds(2)/xwant
            abpolr=-pds(3)/xwant
            kk=(i-2)*4+k
            delxn=-dely
            delyn=delx
            bpoln=(abpolr*delxn+abpolz*delyn)/sqrt(delxn**2+delyn**2)
            btor=fpol(mw)/xwant
            tanbn(kk)=bpoln/btor
c     angbn(kk)=atan2d(ywant,(rangle-xwant))
c     if (angbn(kk).lt.0.0) angbn(kk)=angbn(kk)+360.
            sisibn(kk)=(pds(1)-simag)/(psibry-simag)
c     write (99,*) xwant,ywant,angbn(kk),tanbn(kk),bpoln,btor
               
 1005 continue
            
c     angbn(kk+1)=360.
      tanbn(kk+1)=tanbn(1)
      sisibn(kk+1)=sisibn(1)
      
      continue
c     leaving WEQDSK
      
      return
      
 2030 continue
      print *, ' **G EQDSK file is missing: ', ifnamein
      stop
      
      end
      
