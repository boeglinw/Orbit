      character*80 line
      line = '   hello there   '
      write (6,*) 'no trim:',line
      line = adjustl(line)
      write(6,*) 'trim:',line(1:len_trim(line))
c   
      stop
      end
