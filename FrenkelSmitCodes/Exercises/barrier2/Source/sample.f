      Subroutine Sample(Switch,Ttt)
      Implicit None

      Include 'system.inc'

      Integer          Maxx,I,Ttt,Switch
      Parameter(Maxx=50000)

      Double Precision Kt(Maxx),Ks(Maxx)
      Save Kt,Ks

      If (Switch.Eq.1) Then

C     Initialize Everything

         Do I=1,Maxx
            Kt(I) = 0.0d0
            Ks(I) = 0.0d0
         Enddo

      Elseif(Switch.Eq.2) Then

Ccccccccccccccccccccccccccccccccccccccccc
C     Sample Results: Equation 11.2.1   C
C     Pag. 255 Frenkel/Smit             C
C                                       C
C     The Maximum Time Equals Maxx      C
Ccccccccccccccccccccccccccccccccccccccccc

C     Begin Modification

C     End Modification

      Else

C     Write Results

         Do I=1,Maxx
            If(I.Le.Nstep) Then

               If(Ks(I).Gt.0.5d0) Write(27,*) Dble(I-1),Kt(I)/Ks(I)

            Endif
         Enddo

      Endif

      Return
      End
