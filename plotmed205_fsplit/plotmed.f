      subroutine plotmed (abpolr,abpolz,abt,btotal)
c**********************************************************************
c**                                                                  **
c**     MAIN PROGRAM:  MHD PLOTTING CODE                             **
c**                                                                  **
c**                                                                  **
c**     SUBPROGRAM DESCRIPTION:                                      **
c**          plotme traces out the plasma outermost surface.  This   **
c**          version implements the plotting using DISSPLA.          **
c**                                                                  **
c**     CALLING ARGUMENTS:                                           **
c**                                                                  **
c**     REFERENCES:                                                  **
c**          (1)                                                     **
c**          (2)                                                     **
c**                                                                  **
c**     RECORD OF MODIFICATION:                                      **
c**          12/04/84..........first created                         **
c**          21/11/86..........revised                               **
c**          04/08/96..........revised by Q.Peng                     **
c**          29/01/99..........revised by D. Darrow for FIGOL        **
c**                                                                  **
c**********************************************************************

      include 'dcalc.h'
      include 'eparmdu129.h'
      parameter (nslit=4,m3d=1000)
      parameter (npitch=8*nlimit)
      parameter (nelite=200)

      
      common/mcom/rmaxis,zmaxis
      common/cpitch/ipitch,tanbn(npitch),angbn(npitch),sisibn(npitch)
c     common/cwork3/c(2,nw,nh2),wk(nwrk),copy(nw,nh)
      common/cwork3/c(kubicx,lubicx,kubicy,lubicy),wk(nwrk),
     .     copy(nw,nh),bkx(lubicx+1),bky(lubicy+1),
     .     lkx,lky
      common/csolov/salpha,sbeta,srm,rcentr,bcentr,plasma,dpsimx,ssaa
     .     ,ssee,sbetat,sbetap,sli,svolum,srgeom,strian,isss,jsss
     .     ,risss,zisss
      common/parame/volp(nw),pprime(nw),pres(nw),ffprim(nw),fpol(nw)
      common/consta/pi,tmu,twopi,ioption,ipress,idolim,isurf
      common/gtable/rgrid(nw),zgrid(nh)
      common/cgrid/darea,drgrid,dzgrid,qmaxis,cratio,dfsqe
      common/cpsi/psi(nwnh),psibry,simag,sidif,xpsi(nwnh),eouter
      common/input1/icondn,itek,kdata,ico2,itrace,ierchk,iconvr,ixray
      common/limmm/xlmin,xlmax,ylmin,ylmax,radum
      common/nio/nin,nout,ntty,nrsppc,nrspfc,nttyo,neqdsk,nffile,nsave
      common/limite/limitr,xlim(nlimit),ylim(nlimit)
      common/outp1/mfvers(2)
      common/cxray/rxray(nslit),zxray(nslit),xangle(nangle,nslit)
      common/czero/zero(nwnh),www(nwnh),iweigh
      common/cffol/qpsi(nw),bfpol(nw),cfpol(nw),dfpol(nw),mwfpol,
     .     rbdry(mbdry),zbdry(mbdry),nbdry,xxxsi(nw),sifm,sifb
     .     ,rzero,bzero,piii,qout95,pcurrt(nwnh)
      common/input4/ifname(ntime),islve
      common/ccase/case(6)
      common/dline0/mw2,mh2,idline,siline,br,bz,btor
      common/vtor/presw(nw),preswp(nw),kvtor,rvtor,cwrmid(nw)
     .     ,prw(nw),wpsi(nwnh),presst(nwnh),pressw(nwnh)
c     common/cww/cw(2,nw,nh2),wkw(nwrk),copyw(nw,nh)
      common/cww/cw(kubicx,lubicx,kubicy,lubicy),wkw(nwrk),
     .     copyw(nw,nh),bwx(lubicx+1),bwy(lubicy+1),
     .     lwx,lwy
c     common/cww2/cw2(2,nw,nh2),wkw2(nwrk),copyw2(nw,nh)
      common/cww2/cw2(kubicx,lubicx,kubicy,lubicy),wkw2(nwrk),
     .     copyw2(nw,nh),bwx2(lubicx+1),bwy2(lubicy+1),
     .     lwx2,lwy2
      common/cerror/delerr
      common/mw1com/mw1,mh1
      
      namelist/bdry/nbdry,rbdry,zbdry
      namelist/inrz/rbpol,zbpol,nbpol
      namelist/outbrz/nbpol,rbpol,zbpol,bpolr,bpolz
      character*5 mfvers
      character*35 ifname,filenm,eqdsk,elitefile
      character uday*10,mfitpop*12,case*10
      dimension x3d(m3d),y3d(m3d),z3d(m3d)
      dimension xlims(5),ylims(5),pr(nw),pr2(nw)
      dimension prp(nw),prp2(nw),ffp(nw),ffp2(nw)
      dimension xout(npoint),yout(npoint),zer2(nwnh),
     .     bpolrz(npoint),btorrz(npoint)
      dimension itext(30),workaa(nw),workbb(nw),xiter(kxiter),
     .     itern(kxiter),errort(kxiter),chisq(kxiter)
     .     ,rgri2(nw),zgri2(nh),mfitpopeq(3),pds(6)
      dimension workaw(nw),workbw(nw),workpw(nw),workqw(nw)
      dimension rbpol(npoint),zbpol(npoint),prest2(nwnh)
      dimension bpolr(npoint),bpolz(npoint)
      dimension vars(3),derivs(3),igearwk(3),gearwk(52)
      dimension sia(nelite),sa(nelite),mma(nelite),qa(nelite),
     .     pxa(nelite),
     .     bfppp(nelite),cfppp(nelite),dfppp(nelite)
      dimension sigri1(nw),qpsi1(nw)
      dimension sigri2(nw),qpsi2(nw)
c     external dline,dum1
      data istrpl/0/,nsplot/4/,psitol/1.0e-04/,limid/104/
      data n111/1/,n333/3/
      equivalence (mfitpopeq,mfitpop)

      pi=3.1415926535897932
      imfit1=1
      imfit2=1
      limfag=2
      kbound=0
      m20=0
c
c     store position in xwant and ywant where
c     xwand = r(1) = R
c     ywant = r(3) = Z
c
      xwant = r(1) 
      ywant = r(3) 
d     write (8, *) ' Entered PLOTMED with R= ', (r(i), i=1,3),
d     .  ' xwant=', xwant, ' ywant=', ywant
      
c
c check if position is outside of grid      
      if (xwant.lt.rgrid(1)) go to 65
      if (xwant.gt.rgrid(mw1)) go to 65
      if (ywant.lt.zgrid(1)) go to 65
      if (ywant.gt.zgrid(mh1)) go to 65
c 
c     call dbcevl(rgrid,mw1,zgrid,mh1,c,nw,xwant,ywant,pds,ier)
c     interpolate psi at R,Z and its derivatives
c
      call seva2d(bkx,lkx,bky,lky,c,mw1,mh1,xwant,ywant,pds,ier,n333)

c     calculate the poloidal magnetic field
 
      abpolz=pds(2)/xwant
      abpolr=-pds(3)/xwant
      abpol=sqrt(abpolr*abpolr+abpolz*abpolz)

d     write (8, *) ' abpolz=', abpolz, ' abpolr=', abpolr, ' abpol=', 
d     .  abpol, ' sifb=', sifb, ' sifm=', sifm, ' pds(1)=', pds(1)
      
      
      if (abs(sifb-sifm).gt.1.e-8) then
         xsinow=(pds(1)-sifm)/(sifb-sifm)
      else
         xsinow=1.2
      endif


c     calculate toroidal field component
      
c     Check to see whether sample point is inside limiter region
      call zlim(zeronow,1,1,nbdry,rbdry,zbdry,xwant,ywant,limfag)
      
d     write (8, *) ' xsinow=', xsinow, ' zeronow=', zeronow, ' piii=',
d     .  piii
      
      if (abs(piii).gt.10.) then
         if (zeronow.ge.0.5.and.xsinow.le.1.) then
            
c     take this branch if inside limiter and psi small enough
            
            fpnow=seval(mwfpol,xsinow,xxxsi,fpol,bfpol,cfpol,dfpol)
         else
            fpnow=fpol(mwfpol)
            rhox=999.0
         endif
      else
         fpnow=rzero*bzero
      endif

      
d     write (8, *) ' fpnow=', fpnow
c     print*,'currently in subroutine plotmed'
c     print*,'angle= ',angle
      
      abt=fpnow/xwant

      variable= (angle*3.14)/180.0
      ba=abpolr*cos(variable)+abpolz*sin(variable)
c     print*,'ba= ',ba
      
      
c----------------------------------------------------------------------
c--   Total B field and second harmonic                                --
c--   btotal and bfreq2                                        --
c----------------------------------------------------------------------
      
      
      btotal=sqrt(abt*abt+abpol*abpol)
d     write (8, *) ' Leaving PLOTMED with B= ', 
d     .    abpolr, abt, abpolz, btotal
      
      return
      
 65   continue
      
c     come here if location requested is outside grid
      print *, ' Location requested, (',r(1), r(3), ') is out of grid.'
      
c WB outside grid, use an estimate for the toroidal field

      abpolr = 1e-8
      abpolz = 1e-8
      abt = 1.e-8
      btotal = 2.e-8

      return
      
      end
      
