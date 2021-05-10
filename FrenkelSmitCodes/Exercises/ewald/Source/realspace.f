      Subroutine Realspace(Ureal)
      Implicit None

      Include 'system.inc'

Cccccccccccccccccccccccccccccccccccccccccccc
C     Real Path; Also Direct Calculation   C
C     Loop Over All Particle Pairs         C
Cccccccccccccccccccccccccccccccccccccccccccc

      Integer I,J
      Double Precision Dx,Dy,Dz,Ureal,R,R2,Ir,Dderfc

      Ureal = 0.0d0

C     Start Modification
C     Cameron F Abrams cfa22@drexel.edu 2021 CHET580
      Do I = 1,Npart - 1
         Do J = I + 1,Npart

C     Calculate Distance And Perform Periodic
C     Boundary Conditions

            Dx = Rx(I) - Rx(J)
            Dy = Ry(I) - Ry(J)
            Dz = Rz(I) - Rz(J)

            If (Dx.Gt.Hbox) Then
               Dx = Dx - Box
            Elseif (Dx.Lt. - Hbox) Then
               Dx = Dx + Box
            Endif

            If (Dy.Gt.Hbox) Then
               Dy = Dy - Box
            Elseif (Dy.Lt. - Hbox) Then
               Dy = Dy + Box
            Endif

            If (Dz.Gt.Hbox) Then
               Dz = Dz - Box
            Elseif (Dz.Lt. - Hbox) Then
               Dz = Dz + Box
            Endif

            R2 = Dx*Dx + Dy*Dy + Dz*Dz
            R = sqrt(R2)
            Ir = 1.0/R
            Ureal = Ureal + Z(I)*Z(J)*Ir*Derfc(Alpha*R)
         end do
      end do
C     End   Modification

      Return
      End
