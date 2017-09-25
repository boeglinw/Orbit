
!
! f90 version
!
! can be imported into Python
!

module read_eqdsk_mod
   implicit none
!
   character*10 :: case_t(6) 
   
   integer, parameter :: MB_PTS = 1500
   integer, parameter :: ML_PTS = 1500
   integer, parameter :: neqdsk = 20

   integer :: nw, nh
   integer :: nbbbs,limitr
   
   real :: rdim,zdim,rcentr,rleft,zmid,bcentr
   real :: current,simag,xdum,rmaxis
   real :: zmaxis,sibry   
   
   real,dimension(:,:), allocatable :: psirz
   real,dimension(:), allocatable :: fpol
   real,dimension(:), allocatable :: pres
   real,dimension(:), allocatable :: ffprim
   real,dimension(:), allocatable :: pprime
   real,dimension(:), allocatable :: qpsi
   
   real, dimension(MB_PTS) :: rbbbs
   real, dimension(MB_PTS) :: zbbbs 
   real, dimension(ML_PTS) :: rlim
   real, dimension(ML_PTS) :: zlim
   
contains
      
subroutine read_eqdsk(eq_file)      
   implicit none
   character*80 :: eq_file
   
   integer :: i, j, istat

   integer :: idum
   !
   eq_file = adjustl(eq_file)
   open(neqdsk, file = eq_file(1:len_trim(eq_file)), status = 'old')
   
   read (neqdsk,2000) (case_t(i),i=1,6),idum,nw,nh
   
   ! allocate arrays
   
   allocate( psirz(nw,nh), STAT = istat)
   if (istat .ne. 0) then
      print *, 'Cannot allocate array: psirz ', istat
      return
   endif
   allocate( fpol(nw), STAT = istat)
   if (istat .ne. 0) then
      print *, 'Cannot allocate array: fpol ', istat
      return
   endif
   allocate( pres(nw), STAT = istat)
   if (istat .ne. 0) then
      print *, 'Cannot allocate array: pres ', istat
      return
   endif
   allocate( ffprim(nw), STAT = istat)
   if (istat .ne. 0) then
      print *, 'Cannot allocate array: ffprim ', istat
      return
   endif
   allocate( pprime(nw), STAT = istat)
   if (istat .ne. 0) then
      print *, 'Cannot allocate array: pprime ', istat
      return
   endif
   allocate( qpsi(nw), STAT = istat)
   if (istat .ne. 0) then
      print *, 'Cannot allocate arrays: qpsi ', istat
      return
   endif
   
   read (neqdsk,2020) rdim,zdim,rcentr,rleft,zmid 
   read (neqdsk,2020) rmaxis,zmaxis,simag,sibry,bcentr 
   read (neqdsk,2020) current,xdum,xdum,rmaxis,xdum 
   read (neqdsk,2020) zmaxis,xdum,xdum,xdum,xdum
   read (neqdsk,2020) (fpol(i),i=1,nw) 
   read (neqdsk,2020) (pres(i),i=1,nw) 
   read (neqdsk,2020) (ffprim(i),i=1,nw) 
   read (neqdsk,2020) (pprime(i),i=1,nw) 
   read (neqdsk,2020) ((psirz(i,j),i=1,nw),j=1,nh) 
   read (neqdsk,2020) (qpsi(i),i=1,nw) 
   read (neqdsk,2022) nbbbs,limitr
   read (neqdsk,2020) (rbbbs(i),zbbbs(i),i=1,nbbbs) 
   read (neqdsk,2020) (rlim(i),zlim(i),i=1,limitr)
   ! 
2000 format (6a8,3i4) 
2020 format (5e16.9) 
2022 format (2i5)
   
   close (neqdsk)
   return
   
 end subroutine read_eqdsk
 
 subroutine free_all
   ! deallocate arrays
   deallocate(psirz)
   deallocate(fpol)
   deallocate(pres)
   deallocate(ffprim)
   deallocate(pprime)
   deallocate(qpsi)
   return
 end subroutine free_all
 
end module read_eqdsk_mod
