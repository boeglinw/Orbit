      implicit real*8 (a-h, o-z)
      common /accept_cmn/acc_1, acc_2, acc_3, acc_4, acc_5
      print *, 'enter :  xc, xd, d, xal '
      read (5,*) xc, xd, d, xal 
      faccept = accept(xc, xd, d, xal)
      print *, 'main_program accept = ', faccept
      print *, 'acc_1 = ', acc_1
      print *, 'acc_2 = ', acc_2
      print *, 'acc_3 = ', acc_3
      print *, 'acc_4 = ', acc_4
      print *, 'acc_5 = ', acc_5
      stop
      end
