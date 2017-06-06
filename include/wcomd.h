      common/contor/s1(ntime),s2(ntime),s3(ntime),bpolav(ntime)
      common/hist/eout(ntime),rout(ntime),zout(ntime),doutu(ntime)
     .  ,doutl(ntime),aout(ntime),vout(ntime),betat(ntime),otop(ntime)
     .  ,betap(ntime),ali(ntime),oleft(ntime),oright(ntime),qsta(ntime)
     .  ,rcurrt(ntime),zcurrt(ntime),qout(ntime),olefs(ntime)
     .  ,orighs(ntime),otops(ntime),sibdry(ntime),areao(ntime)
     .  ,wplasm(ntime),elongm(ntime),qqmagx(ntime),terror(ntime)
     .  ,rmagx(ntime),zmagx(ntime),obott(ntime),obots(ntime)
     .  ,alpha(ntime),rttt(ntime),dbpli(ntime),delbp(ntime)
     .  ,rseps(2,ntime),zseps(2,ntime),sepexp(ntime),shearb(ntime)
     .  ,xtch(ntime),ytch(ntime),qpsib(ntime),vertn(ntime),aaq1(ntime)
     .  ,aaq2(ntime),aaq3(ntime),btaxp(ntime),btaxv(ntime)
     .  ,limloc(ntime),simagx(ntime),jerror(ntime),seplim(ntime)
     .  ,wbpol(ntime),taumhd(ntime),betapd(ntime),betatd(ntime)
     .  ,alid(ntime),wplasmd(ntime),taudia(ntime),wbpold(ntime)
     .  ,qmerci(ntime),slantu(ntime),slantl(ntime),zeff(ntime),
     .  zeffr(ntime),tave(ntime),rvsin(ntime),zvsin(ntime),
     .  rvsout(ntime),zvsout(ntime),wpdot(ntime),wbdot(ntime),
     .  vsurfa(ntime),tsaisq(ntime),time(ntime),psiref(ntime),
     .  fluxx(ntime),tavem(ntime),cjor95(ntime),pp95(ntime),
     .  ssep(ntime),yyy2(ntime),xnnc(ntime)
     .  ,cprof(ntime),oring,cjor0(ntime),qqmin(ntime),
     .   chigamt(ntime),ssi95(ntime),rqqmin(ntime)
     .  ,cjor99(ntime),cj1ave(ntime)
      common/cvalue/csilop(nsilop,ntime),crogow(nrogow,ntime),
     .  cmpr2(magpri,ntime),cpasma(ntime),xndnt(ntime)
     . ,cbetap,cli,cqqxis,cbetat,ci0,cipmp2
     . ,ccbrsp(nfcoil,ntime)
      common/exdat2/bcentr(ntime),rcentr,rcencm
      common/exdat3/eccurt(ntime,nesum),pbinj(ntime),vloopt(ntime)
     .   ,pasmat(ntime)
      common/dlc/dfluxc(ntime),sigdlc,rspdlc(nffcur),cdflux(ntime)
      common/comco2/rco2r(nco2r,ntime),rco2v(nco2v,ntime),chordv(nco2v)
     .     ,chordr(nco2r),zcentr,dco2r(ntime,nco2r),dco2v(ntime,nco2v)
      common/outp1/mfvers(2),sbpp,uday,qmflag,header,mco2v,mco2r,
     .             lflag,jflag(ntime),case(6)
      common/parame/volp(nw),pprime(nw),pres(nw),ffprim(nw),fpol(nw)
     .     ,qpsi(nw),r2surf(nw),rpres(nw),psirz(nw,nh),rgrid(nw),
     .      zgrid(nh),ktime,ishot,itime,idum
      common/geqdsk/mw,mh,xdim,zdim,rzero,zmid,rmaxis,zmaxis,ssimag,
     .              ssibry,xlim(nlimit),ylim(nlimit),
     .              rbdry(npoint),zbdry(npoint)
      common/tsrz/zuperts(ntime),zlowerts,rmajts
      character*5 mfvers
      character uday*10,qmflag*3,header*42,case*10
