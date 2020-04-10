      function ifread(comt,test)
      character comt*(*),test*(*),answ*80
      logical ifread
      call nospac(comt,i,j)
10    write (6,'(1x,a,1x,$)') comt(i:j)
      read (5,'(a)') answ
      call nospac(answ,ia1,ia2)
      call nospac(test,it1,it2)
      if (ia2.lt.ia1) goto 10
      ifread=index(test(it1:it2),answ(ia1:ia2)).gt.0
      return
      end
c
      subroutine readin(comt,x)
      character comt*(*)
      call nospac (comt,i,j)
10    write (6,'(1x,a,1x,$)') comt(i:j)
      read (5,*,err=10) x
      return
      end

      subroutine dreadin(comt,x)
      real*8 x
      character comt*(*)
      call nospac (comt,i,j)
10    write (6,'(1x,a,1x,$)') comt(i:j)
      read (5,*,err=10) x
      return
      end
c
      subroutine chread(comt,chr,k,l)
      character comt*(*),chr*(*)
      call nospac(comt,i,j)
10    write (6,'(1x,a,1x,$)') comt(i:j)
      read (5,'(a)') chr
      call nospac(chr,k,l)
      if (l.lt.k) goto 10
      return
      end
c
      subroutine nospac(ch,i1,i2)
c     find first and last non-space characters in an arbitrary length
c     string.  if string is only spaces, i2 will be 1 less than i1
      character ch*(*)
      i2=len(ch)
      i1=1
6000  if (i1.lt.i2.and.ch(i1:i1).eq.' ') then
          i1=i1+1
      goto 6000
      endif
6005  if (i2.ge.i1.and.ch(i2:i2).eq.' ') then
          i2=i2-1
      goto 6005
      endif
      return
      end

      integer function len_trim(ch)

c give the length of the string to the end of the
c last non-black character      

      character*(*) ch

      call nospac(ch, i1,i2)
      if (i1.gt.i2) then
         len_trim = 0
      else
         len_trim = i2
      endif

      return
      end


c
c these are functions for interactive I/O
c
      function rread(comt)
      character comt*(*)
      call nospac (comt,i,j)
10    write (6,'(1x,a,1x,$)') comt(i:j)
      read (5,*,err=10) x
      rread=x
      return
      end

      real*8 function dread(comt)
      real*8 x
      character comt*(*)
      call nospac (comt,i,j)
10    write (6,'(1x,a,1x,$)') comt(i:j)
      read (5,*,err=10) x
      dread=x
      return
      end
c
      function iread(comt)
      character comt*(*)
      call nospac (comt,i,j)
10    write (6,'(1x,a,1x,$)') comt(i:j)
      read (5,*,err=10) x
      iread=nint(x)
      return
      end
c
	function cread(comt)
      character comt*(*),cread*(*), chr*80
      call nospac(comt,i,j)
10    write (6,'(1x,a,1x,$)') comt(i:j)
      read (5,'(a)') chr
      call nospac(chr,k,l)
      if (l.lt.k) goto 10
	cread = chr(k:l)
      return
      end









