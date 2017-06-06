      subroutine write_orbit(effic, accept, Sdl, SdV,  n, jdet)
      include 'dcalc.h'

      include 'eparmdu129.h'
c
c write the coordinates and psi of the current orbit into a file
c
      common /iounit/iorbit, icounter, orbit_fname
      character *80 orbit_fname

c detector index corresponding to the track
      
      integer jdet
      
      parameter (mxorpt=64000)
      common/data5/x(mxorpt+1),y(mxorpt+1),z(mxorpt+1),
     >     x2(mxorpt+1),y2(mxorpt+1), b_r(mxorpt+1), b_phi(mxorpt+1), 
     >     b_z(mxorpt+1)
      
      common /data6/ iflxplt, steppsi(mxorpt+1), pphi(mxorpt+1),
     >     mu(mxorpt+1), energy(mxorpt+1), vpar(mxorpt+1), vperp(mxorpt+1),
     >     bmod(mxorpt+1), rho(mxorpt+1), lborho(mxorpt+1)

c cumulatice sum of sdl, used for checks

      common /sdl_data/sdl_cs(mxorpt+1)

c open filename

      open(iorbit, file = orbit_fname, err = 999)
      
c write header information
      write (iorbit, *) '# ORBIT output '
      write (iorbit, *) '# begin of parameter section '
      write (iorbit, *) '#\ effic = ', effic
      write (iorbit, *) '#\ accept = ', accept
      write (iorbit, *) '#\ stepsize = ', S
      write (iorbit, *) '#\ Sdl = ', Sdl
      write (iorbit, *) '#\ SdV = ', Sdv
      write (iorbit, *) '#\ port_th = ', port
      write (iorbit, *) '#\ port_ph = ', portph
      write (iorbit, *) '#\ detector_index= ', jdet
      write (iorbit, *) '#\ detector_id = ', detector_id(jdet)
      write (iorbit, *) '#\ channel_number = ', channel_number(jdet)
      write (iorbit, *) '# end of parameter section' 
      write (iorbit, *) '#! r[f,0]/ phi[f,1]/ z[f,2]/ x[f,3]/ y[f,4]/ psi[f,5]/ psirel[f,6]/ Em[f,7]/ Sdlcs[f,8]/'//
     >    ' br[f,9]/ bphi[f,10]/ bz[f,11]/'
      do i = 1, n
         j = i+1
         Em_val = Em( get_psirel(steppsi(j)), x(j), z(j) )
         write (iorbit, *) x(j), y(j), z(j), x2(j), y2(j), steppsi(j), get_psirel(steppsi(j)), Em_val, sdl_cs(j), 
     >      b_r(j), b_phi(j), b_z(j)
      enddo
      close(iorbit)
      return
 999  write (6,*) 'cannot open file : ', orbit_fname
      return
      end

     
      
