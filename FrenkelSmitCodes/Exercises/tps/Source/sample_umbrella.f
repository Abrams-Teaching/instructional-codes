      Subroutine Sample_Umbrella(Switch)
      Implicit None

Cccccccccccccccccccccccccccccccccccccccccc
C     Umbrella Sampling Of P(Lambda,T)   C
C     Calculate The Probability That A   C
C     Molecule Is At A Given Position    C
C     Provided That It Is In The Slice   C
Cccccccccccccccccccccccccccccccccccccccccc

      Include 'maxarray.inc'
      Include 'mdtop.inc'
      Include 'mdvelo.inc'
      Include 'traject.inc'

      Logical          Lslice
      Integer          Switch,I
      Double Precision Gg1(Maxumbrel),Order,Gg2,Dddx

      Save Dddx,Gg1

      If(Switch.Eq.1) Then

         Do I=1,Maxumbrel
            Gg1(I) = 0.0d0
         Enddo
      
         Dddx = 1.0d0/Maxdslice
         
      Elseif(Switch.Eq.2) Then

         Do I=1,Natom
            Rxx(I) = Xxold(I,Nslice)
            Ryy(I) = Yyold(I,Nslice)
         Enddo

         Call In_Slice(Lslice,Order)

         If(Lslice) Then
            I = Idint(Dddx*Dsqrt(Order))

            If(I.Lt.1.Or.I.Gt.Maxumbrel) Then
               Write(6,*) 'Order Is Out Of Range !!!'
               Write(6,*) 'Error Sample_Umbrella !!!'
               Write(6,*) 'Order : ',Order
               Write(6,*) 'I     : ',I
               Call Exitt(2)
            Else
               Gg1(I) = Gg1(I) + 1.0d0
            Endif
         Else
            Write(6,*) 'Path Is Not In Slice  !!!'
            Write(6,*) 'Error Sample_Umbrella !!!'
            Call Exitt(2)
         Endif
      
      Else

         Gg2 = 0.0d0

         Do I=1,Maxumbrel
            If(Gg1(I).Gt.0.3d0) Gg2 = Gg2 + Gg1(I)
         Enddo

         Open(32,File='umbrella.dat')
               
         Do I=1,Maxumbrel
            If(Gg1(I).Gt.0.3d0) Then
               Write(32,'(I7,E20.10)') I,Gg1(I)/Gg2
            Endif
         Enddo

         Close(32)
      Endif

      Return
      End
