      Subroutine Ener(Uconf,Itemp)
      Implicit None

      Include 'system.inc'

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Calculate The Total Energy Of A Configuration   C
C     Not Very Efficient, But It Will Do Here         C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Integer I,J,Itemp
      Double Precision Uconf,R2

      Uconf = 0.0d0

      Do I=1,8
         Do J=(I+1),9
            R2 = (Rx(I,Itemp)-Rx(J,Itemp))**2 + 
     &           (Ry(I,Itemp)-Ry(J,Itemp))**2

            If(R2.Lt.1.0d0) 
     &           Uconf = Uconf + (Dsqrt(R2)-1.0d0)**2
         Enddo
      Enddo

      Return
      End
