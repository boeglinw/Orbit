      implicit real*8 (a-h, o-z)

      parameter( n_par = 20)
      real*8 theta_port(n_par), phi_port(n_par), gyro_angle(n_par), pitch_angle(n_par)
      real*8 rdist(n_par),ZDist(n_par), phdangle(n_par)

      real*8, dimension(n_par) :: ports, portsph, al0s, be0s, rds, zds, phda

      integer detec, detector_number(n_par), channel_number(n_par), detector_id(n_par)

c     current position, velocity and magnetic field       

      common/chrcom/r(4),v(4),b(4),s,sstp, port, ports, portph, portsph,
     > efflim, tol, detec, detector_number, channel_number, detector_id

      common/detcom/d,xc,yc,al0,be0,al0s,be0s,rd,rds,phd,phda, xd,yd,zd,zds,nx,ny

