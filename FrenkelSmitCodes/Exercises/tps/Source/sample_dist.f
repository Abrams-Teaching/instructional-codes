      Subroutine Sample_Dist(Switch)
      Implicit None

      Include 'maxarray.inc'

      Logical          Lslice
      Integer          Switch,I
      Double Precision Gg1(Maxumbrel),Gg2,Dddx,D1

      Save Gg1,Dddx

      If(Switch.Eq.1) Then

         Do I=1,Maxumbrel
            Gg1(I) = 0.0d0
         Enddo
               
         Dddx = 1.0d0/Maxdslice
         
      Elseif(Switch.Eq.2) Then

         Call In_Slice(Lslice,D1)

         I = Idint(Dddx*Dsqrt(D1))

         If(I.Lt.1.Or.I.Gt.Maxumbrel) Then
            Write(6,*) 'Order Is Out Of Range !!!'
            Write(6,*) 'Error Sample_Dist !!!'
            Write(6,*) 'I     : ',I
            Call Exitt(2)
         Else
            Gg1(I) = Gg1(I) + 1.0d0
         Endif
            
      Else

         Gg2 = 0.0d0
               
         Do I=1,Maxumbrel
            If(Gg1(I).Gt.0.3d0) Gg2 = Gg2 + Gg1(I)
         Enddo
         
         Open(32,File='distri.dat')

         Do I=1,Maxumbrel
            If(Gg1(I).Gt.0.3d0)
     &           Write(32,'(2E20.10)') Dble(I)*Maxdslice,Gg1(I)/Gg2
         Enddo

         Close(32)
      Endif

      Return
      End
