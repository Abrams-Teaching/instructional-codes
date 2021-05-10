      Subroutine Widom(Weight)
      Implicit None

      Include 'maxarray.inc'
      Include 'system.inc'
      
      Logical Lold

      Double Precision Weight,U,Xf(Maxlength),
     &     Yf(Maxlength),Zf(Maxlength)

Ccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Grow A Test Chain Or Calculatet The Weight    C
C     Of A Chain That Already Exists                C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Lold = Lchain

      Call Grow(Weight,U,Xf,Yf,Zf,Lold)

      Return
      End
