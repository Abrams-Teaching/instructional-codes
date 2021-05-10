      Subroutine Sample_Avhb(Switch)
      Implicit None

Ccccccccccccccccccccccccccccccccccccccccc
C     Sample The Distribution H_B(T)    C
C                                       C
C     Calculate The Probability That A  C
C     Path Is In B At Time t, Provided  C
C     It Is At Least Once In B During   C
C     [0,T] (T>T)                       C
Ccccccccccccccccccccccccccccccccccccccccc

      Include 'maxarray.inc'
      Include 'mdtop.inc'
      Include 'traject.inc'

      Integer          Switch,J
      Logical          Lready
      Double Precision Gg1(Maxtraject),Gg2
      
      Save Gg1,Gg2

      If(Switch.Eq.1) Then

         Gg2 = 0.0d0

         Do J=1,Maxtraject
            Gg1(J) = 0.0d0
         Enddo
      
      Elseif(Switch.Eq.2) Then
      
         Gg2 = Gg2 + 1.0d0
         
         Lready = .False.

         Do J=1,Nslice
            If(Lb_Old(J)) Then
               Gg1(J) = Gg1(J) + 1.0d0
               Lready = .True.
            Endif
         Enddo
         
         If(.Not.La_Old(1)) Then
            Write(6,*) 'Error ! Initial Path Is Not In A !!!'
            Write(6,*) 'Error In Sample_Avhb !!!'
            Call Exitt(2)
         Elseif(.Not.Lready) Then
            Write(6,*) 'Error ! Path Is Never In B !!!'
            Write(6,*) 'Error In Sample_Avhb !!!'
            Call Exitt(2)
         Endif
         
      Elseif(Switch.eq.3) Then

         Open(32,File='avhb.dat')

         Do J=1,Nslice
            If(Gg2.Gt.0.5d0) Then
               Write(32,'(2E20.10)')  Tstep*Dble(Nshort*J),Gg1(J)/Gg2
            Else
               Write(32,'(E20.10,A)') Tstep*Dble(Nshort*J),'   0.0d0'
            Endif
         Enddo

         Close(32)
      Endif

      Return
      End
