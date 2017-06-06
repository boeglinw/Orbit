c function to calculate the accpentance of a detectro coll. configuration
c in one dimension
      real*8 function accept(theta)
      implicit real*8 (a-h, o-z)
c      common /accept_cmn/acc_x1, acc_x2, acc_x3, acc_x4, acc_x5
c
c INPUT: rd : detector radius
c         d : distance coll. - detector      
c

	theta_max = atan(2*rd/d)

	cc = d*tan(theta)
      if (xc .le. xd) then
         full_opening = 2.*xc
         t_th1p = (xd - xc - xal)/d
         t_th1n = (xc - xd - xal)/d
      else
         full_opening = 2.*xd
         t_th1n = (xd - xc - xal)/d
         t_th1p = (xc - xd - xal)/d
      endif
      t_th2p = (xc + xd - xal)/d
      t_th2n = -(xc + xd + xal)/d
      th1p = atan(t_th1p)
      th1n = atan(t_th1n)
      th2p = atan(t_th2p)
      th2n = atan(t_th2n)
      acc_x1 =     full_opening*(sin(th1p) - sin(th1n))
      acc_x2 =     (xc + xd + xal)*(sin(th1n) - sin(th2n))
      acc_x3 =  -d*(cos(th1n) - cos(th2n)) 
      acc_x4 =   (xc + xd - xal)*(sin(th2p) - sin(th1p)) 
      acc_x5 =   d*(cos(th2p) - cos(th1p))
      accept =    acc_x1 + acc_x2 + acc_x3 + acc_x4 + acc_x5
      return
      end
