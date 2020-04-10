

      function srchz (z, zlist, nzl)
c*************************************************************
c    This function searches through a list of z values (assumed
c  sorted lowest to highest) to find the pair of z values in the list
c  between which the given z lies.  The routine returns the index of
c  the first of the pair of z's.  If the z value is out of range of
c  the list, then an index of -1 is returned.
c*************************************************************

        implicit real*8 (a-h, o-z)

        save
	integer nzl, top, bot, mid
        integer srchz
	real*8 z, zlist(nzl)

d       write(8,*) ' SRCHZ entered, z= ', z, ' zlist(1),(n)=',
d    .    zlist(1), zlist(nzl)

        if ((z .gt. zlist(nzl)) .or. (z .lt. zlist(1))) then
          srchz = -1
          return
        endif

d       write (8, *) ' In SRCHZ, passed out-of-bounds test'

c   set up binary search

	top = nzl
        bot = 1
        mid = (nzl + 1) / 2

c    WHILE loop follows

        do while ((mid .ne. top) .and. (mid .ne. bot))

d         write (8, *) ' top=', top, ' mid=', mid, ' bot=', bot

          if ((z .ge. zlist(mid)) .and. (z .le. zlist(mid+1))) then

c         found the right interval, so return

            srchz = mid

d           write(8,*) ' SRCHZ: satisfied condition inside while loop:',
d    .        srchz, zlist(mid), zlist(mid+1)

            return
          endif
          if (z .lt. zlist(mid)) then
            top = mid
          else             
            bot = mid
          endif
          mid = (top + bot) / 2
        enddo
        srchz = mid
        if (mid .ge. nzl) srchz = nzl - 1

d       write (8,*)' leaving SRCHZ, srchz=', srchz, zlist(srchz),
d    .    zlist(srchz+1)

        return
        end
