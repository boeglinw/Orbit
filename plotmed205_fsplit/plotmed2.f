
      function plotmed2(flx,w,h) 
c**********************************************************************
c**                                                                  **
c**     MAIN PROGRAM:  MHD PLOTTING CODE                             **
c**                                                                  **
c**                                                                  **
c**     SUBPROGRAM DESCRIPTION:                                      **
c**          This is a modification to the PLOTMED program that      ** 
c**          returns the value of the poloidal flux at (w,h)         **
c**                                                                  **
c**     CALLING ARGUMENTS:                                           **
c**          FLX: the poloidal flux at the location (returned)       **
c**          W:   Major radius of point where flux desired (input)   **
c**          H:   Vertical coord of point where flux desired (input) **
c**                                                                  **
c**     RECORD OF MODIFICATION:                                      **
c**          12/04/84..........first created                         **
c**          21/11/86..........revised                               **
c**          04/08/96..........revised by Q.Peng                     **
c**          29/01/99..........revised by D. Darrow for EIGOL        **
c**          05/11/10..........changed to a function by. W Boeglin   **
c**                                                                  **
c**                                                                  **
c**  returns 1 if everything worked fine                             **
c**          0 if outside the grid                                   **
c**                                                                  **
c**                                                                  **
c**********************************************************************

      implicit real*8 (a-h, o-z)

      integer plotmed2

      include 'eparmdu129.h'
      parameter (nslit=4,m3d=1000)
      parameter (npitch=8*nlimit)
      parameter (nelite=200)
      common/mcom/rmaxis,zmaxis
      common/cpitch/ipitch,tanbn(npitch),angbn(npitch),sisibn(npitch)
c     common/cwork3/c(2,nw,nh2),wk(nwrk),copy(nw,nh)
      common/cwork3/c(kubicx,lubicx,kubicy,lubicy),wk(nwrk),
     .         copy(nw,nh),bkx(lubicx+1),bky(lubicy+1),
     .         lkx,lky
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
     .             rbdry(mbdry),zbdry(mbdry),nbdry,xxxsi(nw),sifm,sifb
     .            ,rzero,bzero,piii,qout95,pcurrt(nwnh)
      common/input4/ifname(ntime),islve
      common/ccase/case(6)
      common/dline0/mw2,mh2,idline,siline,br,bz,btor
      common/vtor/presw(nw),preswp(nw),kvtor,rvtor,cwrmid(nw)
     .            ,prw(nw),wpsi(nwnh),presst(nwnh),pressw(nwnh)
c     common/cww/cw(2,nw,nh2),wkw(nwrk),copyw(nw,nh)
      common/cww/cw(kubicx,lubicx,kubicy,lubicy),wkw(nwrk),
     .         copyw(nw,nh),bwx(lubicx+1),bwy(lubicy+1),
     .         lwx,lwy
c     common/cww2/cw2(2,nw,nh2),wkw2(nwrk),copyw2(nw,nh)
      common/cww2/cw2(kubicx,lubicx,kubicy,lubicy),wkw2(nwrk),
     .         copyw2(nw,nh),bwx2(lubicx+1),bwy2(lubicy+1),
     .         lwx2,lwy2
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
     .          bpolrz(npoint),btorrz(npoint)
      dimension itext(30),workaa(nw),workbb(nw),xiter(kxiter),
     .  itern(kxiter),errort(kxiter),chisq(kxiter)
     .  ,rgri2(nw),zgri2(nh),mfitpopeq(3),pds(6)
      dimension workaw(nw),workbw(nw),workpw(nw),workqw(nw)
      dimension rbpol(npoint),zbpol(npoint),prest2(nwnh)
      dimension bpolr(npoint),bpolz(npoint)
      dimension vars(3),derivs(3),igearwk(3),gearwk(52)
      dimension sia(nelite),sa(nelite),mma(nelite),qa(nelite),
     .          pxa(nelite),
     .          bfppp(nelite),cfppp(nelite),dfppp(nelite)
      dimension sigri1(nw),qpsi1(nw)
      dimension sigri2(nw),qpsi2(nw)
c      external dline,dum1
      data istrpl/0/,nsplot/4/,psitol/1.0e-04/,limid/104/
      data n111/1/,n333/3/
      equivalence (mfitpopeq,mfitpop)

c----------------------------Executable code starts here-------------------
      pi=3.1415926535897932
      xwant = w 
      ywant = h 
d     write (8, *) ' Entered PLOTMED2 with xwant=', xwant, ' ywant=', ywant


      if (xwant.lt.rgrid(1)) go to 65
      if (xwant.gt.rgrid(mw1)) go to 65
      if (ywant.lt.zgrid(1)) go to 65
      if (ywant.gt.zgrid(mh1)) go to 65
c        call dbcevl(rgrid,mw1,zgrid,mh1,c,nw,xwant,ywant,pds,ier)
      call seva2d(bkx,lkx,bky,lky,c,mw1,mh1,xwant,ywant,pds,ier,n333)
      flx = pds(1)

c end change
d     write (8, *) ' Leaving PLOTMED2 with flx=', flx
      plotmed2 = 1
      return

 65   continue

c     come here if point is out of bounds
      flx = 0.0
      print *, ' In PLOTMED2, (R,Z) is outside grid: (',w,h,')'
d     write (8,*) ' Leaving PLOTMED2, (R,Z) is outside grid, so flx=0'

      plotmed2 = 0
      return

      end
