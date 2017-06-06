c function to calculate the accpentance of a detectro coll. configuration
c in one dimension
      real*8 function trans(rd, rc, d)
      implicit real*8 (a-h, o-z)
      real,parameter:: pi = 3.14159265359 
      common /accept_cmn/acc_x1, acc_x2, acc_x3, acc_x4, acc_x5
c
c INPUT: rd : active area of detector radius
c	 rc : collimator radius	
c         d : distance coll. - detector    
c      theta: angle of entrance  c

	theta_max = atan(2*rc,d)

	cc = d*tan(theta)
	x = cc/2.
	phi = acos(x/rd)
	psi = acos(x/rd)
	y = rd*sin(phi)
	area1 = phi*pi*(rd**2.) - y*x
	area2 = psi*pi*(rd**2) - y*x
	area = area1 + area2
	return
	end
