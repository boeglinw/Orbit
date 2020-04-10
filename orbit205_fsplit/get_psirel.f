      real*8 function get_psirel(psi)

      implicit real*8 (a-h, o-z)

c
c INPUT psi
c
c OUTPUT psirel : relative psi (max. in the center)
c
c           psi: absolute value of psi      

c     /CFFOL/ contains info from EFIT magnetics
      
      include 'eparmdu129.h'
      
      
      common/cffol/qpsi(nw),bfpol(nw),cfpol(nw),dfpol(nw),mwfpol,
     >     rbdry(mbdry),zbdry(mbdry),nbdry,xxxsi(nw),sifm,sifb
     >     ,rzero,bzero,piii,qout95,pcurrt(nwnh)
      
      get_psirel = (psi - sifb) / (sifm - sifb)

      return
      end
