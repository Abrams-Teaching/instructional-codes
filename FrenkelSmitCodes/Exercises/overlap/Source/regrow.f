      Subroutine Regrow(Dispa,Dispb)
      Implicit None

      Include 'maxarray.inc'
      Include 'system.inc'

Ccccccccccccccccccccccccccccccccccc
C     Regrow A Chain Using Cbmc   C
Ccccccccccccccccccccccccccccccccccc

      Integer I
      Double Precision Dispa,Dispb,Unew,Uold,
     &       Weinew,Weiold,Xf(Maxlength),
     &       Yf(Maxlength),Zf(Maxlength),
     &       Ran_Uniform

      Dispb = Dispb + 1.0d0

      Call Grow(Weiold,Uold,Xf,Yf,Zf,.True.)
      Call Grow(Weinew,Unew,Xf,Yf,Zf,.False.)

Ccccccccccccccccccccccccccccc
C     Test For Acceptance   C
Ccccccccccccccccccccccccccccc

      If(Ran_Uniform().Lt.(Weinew/Weiold)) Then

Ccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Accepted !! Update Energy And Coordinates   C
Ccccccccccccccccccccccccccccccccccccccccccccccccccc

         Do I=1,Nlength
            X(I) = Xf(I)
            Y(I) = Yf(I)
            Z(I) = Zf(I)
         Enddo

         Dispa = Dispa + 1.0d0
         Usim  = Usim  + Unew - Uold
      Endif

      Return
      End
